import collections
import logging
import itertools
import os
import pickle
import re
import json
import tempfile
from pathlib import Path
from time import time, perf_counter, process_time
import multiprocessing as mp

import numpy
import tables
from tqdm import tqdm

from ..db import Database
from ..exceptions import DBConsistencyError
from ..models import ProteinEntry
from ..tablefmt import ProteinCacheInfo, RootHOGMetaTable

logger = logging.getLogger(__name__)
Protein = collections.namedtuple("Protein", ("entry_nr", "hog_id", "group"))


def are_orthologous(a: Protein, b: Protein):
    if a.entry_nr == b.entry_nr:
        return False
    for (num1, char1), (num2, char2) in zip(a.hog_id, b.hog_id):
        if num1 != num2:
            return True
        if char1 != char2:
            return False
    return True


def create_job_files(db_path, out_prefix):
    fam_sizes = collections.Counter()
    singletons = []
    re_fam = re.compile(rb"HOG:[A-Z]?(\d+)(.[\w.])?")
    with tables.open_file(db_path) as h5:
        etab: tables.Table = h5.get_node("/Protein/Entries")
        entry_to_fam = numpy.zeros(etab.nrows + 1, dtype=numpy.int32)
        for row in etab:
            if len(row["OmaHOG"]) > 0:
                fam = int(re_fam.match(row["OmaHOG"]).group(1))
                entry_to_fam[row["EntryNr"]] = fam
                fam_sizes.update((fam,))
            else:
                singletons.append(
                    (
                        int(row["EntryNr"]),
                        int(row["OmaGroup"]),
                    )
                )
    with open(out_prefix + "_singleton.pkl", "wb") as fh:
        pickle.dump(["process_singletons", singletons], fh)
    # save entry_to_fam to allow for fast lookup later on
    numpy.save(out_prefix + "_entry_to_fam.npy", entry_to_fam)

    def yield_buckets(counts, nr_elem=10000):
        cur_bucket, cur_size = [], 0
        target = nr_elem**2
        for fam, size in sorted(counts.items(), key=lambda x: -x[1]):
            if size**2 > target:
                for rng in range(0, size, nr_elem):
                    yield [(fam, (rng, min(rng + nr_elem, size)))]
            elif cur_size + size**2 < target and len(cur_bucket) < 2 * nr_elem:
                cur_bucket.append((fam,))
                cur_size += size**2
            else:
                if cur_bucket:
                    yield cur_bucket
                cur_bucket = [(fam,)]
                cur_size = size**2
        if cur_bucket:
            yield cur_bucket

    for job, bucket in enumerate(yield_buckets(fam_sizes)):
        with open(f"{out_prefix}_fam-{job:03d}.pkl", "wb") as fh:
            pickle.dump(["process_family", bucket], fh)
    logger.info("wrote %d family job-files and 1 singleton file", job)


# -----------------------------------------------------
# parallel worker code to generate a temporary hdf5 file with
# all the VPairs sorted by RootHOG to efficiently scan the
# orthologs per Family when processing the jobs for the cachebuilder
#
# Overview of this code:
# ┌───────────────────────────┐
# │ Main process              │
# │  - scans genomes          │
# │  - distributes to workers │
# │  - collects produced chunk paths
# └───────────────────────────┘
#          │
#          ▼
# ┌───────────────────────────┐
# │ Worker process (nproc)    │
# │  - open db                │
# │  - read VPairs tables     │
# │  - write chunks (sorted)  │
# └───────────────────────────┘
#          │
#          ▼
# ┌───────────────────────────┐
# │ GlobalMerger              │
# │  - family-wise merge      │
# │  - global FamilyIndex     │
# └───────────────────────────┘

# --------------------------------------------------------------------
# Buffered chunk writer
# --------------------------------------------------------------------
FILTERS = tables.Filters(complevel=5, complib="blosc2")
DTYPE_VPAIRS = numpy.dtype([("Fam", numpy.int32), ("EntryNr1", numpy.int32), ("EntryNr2", numpy.int32)])
DTYPE_FAMILY_INDEX = numpy.dtype([("Fam", numpy.int32), ("Start", numpy.int64), ("End", numpy.int64)])


class BufferedChunkWriter:
    """Accumulate VPairs rows across multiple genomes, flush to .npy
    only when buffer exceeds target_rows."""

    def __init__(self, out_dir, target_rows=5_000_000_000):
        self.out_dir = Path(out_dir)
        self.target_rows = target_rows
        self.buffer = []
        self.total_rows = 0
        self.chunk_id = 0
        self.out_dir.mkdir(exist_ok=True)

    def add_rows(self, arr):
        self.buffer.append(arr)
        self.total_rows += len(arr)
        if self.total_rows >= self.target_rows:
            return self.flush()
        return None

    def flush(self):
        if not self.buffer:
            return None
        chunk = numpy.concatenate(self.buffer)
        chunk.sort(order=["Fam", "EntryNr1"], kind="mergesort")

        fams, start, count = numpy.unique(chunk["Fam"], return_index=True, return_counts=True)
        fam_index = numpy.zeros(len(fams), dtype=DTYPE_FAMILY_INDEX)
        fam_index["Fam"], fam_index["Start"], fam_index["End"] = fams, start, start + count

        logger.debug("flushing chunk %d with %d rows to %s", self.chunk_id, len(chunk), self.out_dir)
        out_path = self.out_dir / f"chunk_{self.chunk_id:05d}.h5"
        with tables.open_file(out_path, "w", filters=FILTERS) as h5:
            h5.create_table("/", "AllVPairs", obj=chunk)
            h5.create_table("/", "FamilyIndex", obj=fam_index)

        self.buffer.clear()
        self.total_rows = 0
        self.chunk_id += 1
        return str(out_path)


# --------------------------------------------------------------------
# 1. GenomeExtractor (multiprocessing) worker
# --------------------------------------------------------------------
def extract_genome_vpairs(args: tuple[str, list[str], str, int, str]) -> list[os.PathLike]:
    db_h5_path, genomes, out_dir, chunk_size, entry_to_fam_path = args
    writer = BufferedChunkWriter(out_dir, target_rows=chunk_size)
    # Load entry → fam mapping as memory-map
    entry_to_fam = numpy.load(entry_to_fam_path, mmap_mode="r")
    chunk_paths = []

    with tables.open_file(db_h5_path, "r") as h5:
        for genome in genomes:
            tab = h5.get_node(f"/PairwiseRelation/{genome}/VPairs")
            nrows = tab.nrows
            batch_size = 5_000_000
            for start in range(0, nrows, batch_size):
                chunk = tab[start : start + batch_size]
                if len(chunk) == 0:
                    continue
                fams = entry_to_fam[chunk["EntryNr1"]]
                fam_arr = numpy.empty(len(chunk), dtype=DTYPE_VPAIRS)
                fam_arr["Fam"], fam_arr["EntryNr1"], fam_arr["EntryNr2"] = fams, chunk["EntryNr1"], chunk["EntryNr2"]

                path = writer.add_rows(fam_arr)
                if path is not None:
                    chunk_paths.append(path)
    final = writer.flush()
    if final is not None:
        chunk_paths.append(final)
    return chunk_paths


# --------------------------------------------------------------------
# Global merge of all sorted chunks
# --------------------------------------------------------------------
class GlobalMerger:
    def __init__(self, buf_write=100_000_000):
        self.buf_write = buf_write

    def merge_chunks(self, chunk_files, out_h5):
        n_chunks = len(chunk_files)
        fins = [tables.open_file(f, "r") for f in chunk_files]
        tabs = [f.get_node("/AllVPairs") for f in fins]
        family_indices = [f.get_node("/FamilyIndex").read() for f in fins]
        positions = [0] * n_chunks
        tot_vpairs = sum(t.nrows for t in tabs)
        logger.info("Merging %d chunks into %s. %d relevant VPairs in total", n_chunks, out_h5, tot_vpairs)

        with tables.open_file(out_h5, "w", filters=FILTERS) as fout:
            out_tab = fout.create_table("/", "AllVPairs", description=DTYPE_VPAIRS, expectedrows=tot_vpairs)
            fam_index = []
            current_offset = 0
            buffer = []
            total_fams = sum(len(idx) for idx in family_indices)
            pbar = tqdm(total=total_fams, desc="Merging")

            while True:
                next_fam = None
                for ci in range(n_chunks):
                    if positions[ci] < len(family_indices[ci]):
                        fam = family_indices[ci]["Fam"][positions[ci]]
                        if next_fam is None or fam < next_fam:
                            next_fam = fam
                if next_fam is None:
                    break

                # Collect all rows for this family from all chunks
                fam_rows = []
                for ci in range(n_chunks):
                    if positions[ci] < len(family_indices[ci]) and family_indices[ci]["Fam"][positions[ci]] == next_fam:
                        s, e = family_indices[ci]["Start"][positions[ci]], family_indices[ci]["End"][positions[ci]]
                        fam_rows.append(tabs[ci][s:e])
                        positions[ci] += 1

                fam_block = numpy.concatenate(fam_rows)
                fam_block.sort(order="EntryNr1", kind="mergesort")
                buffer.append(fam_block)
                fam_index.append((next_fam, current_offset, current_offset + len(fam_block)))
                current_offset += len(fam_block)

                if sum(len(b) for b in buffer) >= self.buf_write:
                    out_tab.append(numpy.concatenate(buffer))
                    buffer.clear()
                    out_tab.flush()
                pbar.update(1)

            if buffer:
                out_tab.append(numpy.concatenate(buffer))
                out_tab.flush()

            fout.create_table("/", "FamilyIndex", obj=numpy.array(fam_index, dtype=DTYPE_FAMILY_INDEX))
            pbar.close()

        for f in fins:
            f.close()


def build_allvpairs_hdf5(
    db_path: os.PathLike, entry_to_fam_path: os.PathLike, out_path: os.PathLike, nproc=8, target_rows=5_000_000_000
):
    """Build the AllVPairs table from the database."""
    with tables.open_file(db_path, "r") as h5:
        genomes: numpy.ndarray = h5.get_node("/Genome").read(field="UniProtSpeciesCode")
    genomes = numpy.char.decode(genomes, "utf-8")
    numpy.random.shuffle(genomes)
    split_genomes = numpy.array_split(genomes, nproc)
    tmp_dirs = [tempfile.mkdtemp() for _ in range(nproc)]

    args = [(db_path, list(split_genomes[i]), tmp_dirs[i], target_rows, entry_to_fam_path) for i in range(nproc)]
    logger.info("Extracting VPairs in parallel...")
    # Parallel ectraction of VPairs
    with mp.Pool(nproc) as pool:
        chunk_lists = pool.map(extract_genome_vpairs, args)

    raw_chunks = [p for sublist in chunk_lists for p in sublist]
    logger.info(f"{len(raw_chunks)} chunks written.")

    # Merge all chunks into a single table
    merger = GlobalMerger()
    merger.merge_chunks(raw_chunks, out_path)


# END of parallel code for temporary hdf5 file
# ----------------------------------------


def process_job_file(job_file: os.PathLike, db_fpath: os.PathLike, vp_fpath: os.PathLike, out: os.PathLike):
    with open(job_file, "rb") as fh:
        jobdata = pickle.load(fh)
    job, payload = jobdata
    with CacheBuilder(db_fpath, vp_db_path=vp_fpath, out_path=out) as builder:
        func = getattr(builder, job)
        if job == "process_singletons":
            func(payload)
        else:
            for args in payload:
                func(*args)


def log_timing(func):
    def wrapper(*args, **kwargs):
        start_wall = perf_counter()
        start_cpu = process_time()
        result = func(*args, **kwargs)
        end_wall = perf_counter()
        end_cpu = process_time()
        wall = end_wall - start_wall
        cpu = end_cpu - start_cpu
        efficiency = cpu / wall if wall > 0 else 0
        logger.info(f"{func.__name__}: CPU={cpu:.6f}s, Wall={wall:.6f}s, Efficiency={efficiency:.2%}")
        return result

    return wrapper


class CacheBuilder:
    def __init__(self, db_fpath, vp_db_path, out_path):
        self.db_fpath = db_fpath
        self.vp_db_path = vp_db_path
        self.out_path = out_path
        self.db = None
        self.vp = None
        self.h5 = None
        self.cnts = []
        self.json_buffer = []
        self.json_offsets = []
        self._buffer_offset = 0
        self._offset_dtype = [("Fam", "i4"), ("offset", "i8"), ("length", "i4")]

    def __enter__(self):
        self.db = Database(self.db_fpath)
        self.h5 = self.db.get_hdf5_handle()
        self.vp = tables.open_file(self.vp_db_path, "r")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.db.close()
        self.vp.close()
        self.save()

    def process_family(self, fam, rng=None):
        t0 = time()
        members = self.load_fam_members(fam)
        self.cnts.append(self.analyse_fam(fam, members, rng))
        if rng is None or rng[0] == 0:
            fam_json = self.compute_familydata_json(fam, members)
            self.store_familydata_json_result(fam, fam_json)
        logger.info(f"processed family {fam} with {len(members)}, took {time() - t0}sec")

    def process_singletons(self, singletons):
        t0 = time()
        self.cnts.append(self.analyse_singleton(singletons))
        logger.info(f"processed %d singletons in %.2f sec", len(singletons), time() - t0)

    def load_fam_members(self, fam):
        members = []
        vals = {k: self.db.format_hogid(x).encode("utf-8") for k, x in zip(("fam", "fam_next"), (fam, fam + 1))}
        for row in self.h5.get_node("/Protein/Entries").where(
            "(fam <= OmaHOG) & (OmaHOG < fam_next)",
            condvars=vals,
        ):
            subhog_tokens = row["OmaHOG"].decode().split(".")[1:]
            parsed_parts = []
            for part in subhog_tokens:
                match = re.match(r"(\d+)([a-z]+)", part)
                if match:
                    num = int(match.group(1))
                    char = match.group(2)
                    parsed_parts.append((num, char))
                else:
                    raise ValueError(f"Invalid part format in '{part}'")
            members.append(Protein(row["EntryNr"], parsed_parts, row["OmaGroup"]))
        return members

    def load_vps(self, entry_nr):
        return self.db.get_vpairs(entry_nr)["EntryNr2"]

    @log_timing
    def load_vps_for_family(self, fam):
        fam_idx_tab = self.vp.get_node("/FamilyIndex")
        fams = fam_idx_tab.read(field="Fam")
        fam_pos = numpy.searchsorted(fams, fam)
        if fam_pos >= len(fams) or fams[fam_pos] != fam:
            return numpy.zeros((0,), dtype=DTYPE_VPAIRS)
        start, end = fam_idx_tab[fam_pos]["Start"], fam_idx_tab[fam_pos]["End"]
        vps_tab = self.vp.get_node("/AllVPairs")
        return vps_tab[start:end]

    def load_grp_members(self, group):
        return [row["EntryNr"] for row in self.h5.get_node("/Protein/Entries").where(f"OmaGroup == {group}")]

    def analyse_fam(self, fam, fam_members, rng=None):
        logger.debug(f"analysing orthology of family {fam} with {len(fam_members)} members; doing range {rng}")
        grp_members = {grp: set(self.load_grp_members(grp)) for grp in set(z.group for z in fam_members if z.group > 0)}
        nr_memb = len(fam_members)
        fam_iter = iter(fam_members)
        if rng is not None:
            nr_memb = rng[1] - rng[0]
            fam_iter = itertools.islice(fam_members, rng[0], rng[1])

        fam_vps = self.load_vps_for_family(fam)
        counts = numpy.zeros(nr_memb, dtype=tables.dtype_from_descr(ProteinCacheInfo))
        time_vps_wall, time_vps_cpu, time_ind_wall, time_ind_cpu, cpu_0 = 0, 0, 0, 0, process_time()
        for i, p1 in tqdm(
            enumerate(fam_iter),
            disable=len(fam_members) < 500,
            desc=f"Processing family {fam}",
            total=nr_memb,
        ):
            t0_cpu, t0_wall = process_time(), perf_counter()
            vps = set(fam_vps[fam_vps["EntryNr1"] == p1.entry_nr]["EntryNr2"])
            t1_cpu, t1_wall = process_time(), perf_counter()
            ind_orth = set(p2.entry_nr for p2 in fam_members if are_orthologous(p1, p2))
            t2_cpu, t2_wall = process_time(), perf_counter()
            time_vps_wall += t1_wall - t0_wall
            time_vps_cpu += t1_cpu - t0_cpu
            time_ind_wall += t2_wall - t1_wall
            time_ind_cpu += t2_cpu - t1_cpu
            grp = grp_members.get(p1.group, set([])) - {p1.entry_nr}
            counts[i]["EntryNr"] = p1.entry_nr
            counts[i]["NrPairwiseOrthologs"] = len(vps)
            counts[i]["NrHogInducedPWOrthologs"] = len(ind_orth)
            counts[i]["NrHogInducedPWParalogs"] = len(fam_members) - len(ind_orth) - 1
            counts[i]["NrOMAGroupOrthologs"] = len(grp)
            counts[i]["NrAnyOrthologs"] = len(vps | ind_orth | grp)
        cpu_1 = process_time()
        logger.debug(
            "timings for family %s: vps=[wall: %.3f; cpu: %.3f] ind=[wall: %.3f; cpu: %.3f]",
            fam,
            time_vps_wall,
            time_vps_cpu,
            time_ind_wall,
            time_ind_cpu,
        )
        logger.debug("  total cpu time for family %s: %.3f sec", fam, cpu_1 - cpu_0)
        logger.debug(
            "  efficiency vps: %.2f%%, ind: %.2f%%",
            (time_vps_cpu / time_vps_wall * 100) if time_vps_wall > 0 else 0,
            (time_ind_cpu / time_ind_wall * 100) if time_ind_wall > 0 else 0,
        )
        logger.debug(
            "  overall efficiency: %.2f%%",
            ((cpu_1 - cpu_0) / (time_vps_wall + time_ind_wall) * 100) if (time_vps_wall + time_ind_wall) > 0 else 0,
        )
        return counts

    @log_timing
    def analyse_singleton(self, singletons):
        logger.info("analysing %d singletons", len(singletons))
        counts = numpy.zeros(len(singletons), dtype=tables.dtype_from_descr(ProteinCacheInfo))
        for i, (entry_nr, group_nr) in enumerate(singletons):
            vps = set(self.load_vps(entry_nr))
            grp_members = set([])
            if group_nr > 0:
                grp_members = set(self.load_grp_members(group_nr))
            counts[i] = (entry_nr, len(vps), 0, 0, len(grp_members), len(vps | grp_members))
        return counts

    @log_timing
    def compute_familydata_json(self, fam, fam_members):
        famhog_id = self.db.format_hogid(fam)
        logger.debug("family data for %s with %d members", fam, len(fam_members))
        # TODO: enable this again once it's clear if still needed and how to access GO at the time of cache building
        if len(fam_members) > 0:
            # this will likely fail to compute the MDS for so many points
            # let's skip it for now.
            genes_null_similarity = set(p.entry_nr for p in fam_members)
        else:
            try:
                (
                    genes_null_similarity,
                    gene_similarity_vals,
                ) = self.db.get_gene_similarities_hog(famhog_id)
            except Exception as e:
                logger.error("gene_similarity failed for %s: %s", fam, e)
                raise
        final_json_output = []
        for p1 in fam_members:
            to_append = {}
            protein = ProteinEntry(self.db, p1.entry_nr)
            to_append["id"] = p1.entry_nr
            to_append["protid"] = protein.omaid
            to_append["sequence_length"] = protein.sequence_length
            to_append["taxon"] = protein.genome.species_and_strain_as_dict
            to_append["xrefid"] = protein.canonicalid
            to_append["gc_content"] = protein.gc_content
            to_append["nr_exons"] = protein.nr_exons
            if p1.entry_nr in genes_null_similarity:
                to_append["gene_similarity"] = None
            else:
                to_append["gene_similarity"] = gene_similarity_vals[p1.entry_nr]
            final_json_output.append(to_append)
        return json.dumps(final_json_output)

    def store_familydata_json_result(self, fam, result):
        encoded_json = result.encode("utf-8")
        json_as_np = numpy.ndarray((len(encoded_json),), buffer=encoded_json, dtype=tables.StringAtom(1))
        self.json_buffer.append(json_as_np)
        self.json_offsets.append(numpy.array([(fam, self._buffer_offset, len(encoded_json))], dtype=self._offset_dtype))
        self._buffer_offset += len(encoded_json)

    def save(self):
        logger.info(f"writing results to {self.out_path}")
        with tables.open_file(self.out_path, "w") as h5:
            buf = h5.create_earray(
                "/family_json",
                "buffer",
                tables.StringAtom(1),
                (0,),
                createparents=True,
                expectedrows=1e9,
            )
            for el in self.json_buffer:
                buf.append(el)
            buf.flush()
            if len(self.json_offsets) > 0:
                off = numpy.concatenate(self.json_offsets)
            else:
                off = numpy.zeros(0, dtype=self._offset_dtype)
            h5.create_table("/family_json", "offset", obj=off)

            if len(self.cnts) > 0:
                cnts = numpy.concatenate(self.cnts)
            else:
                cnts = numpy.zeros(0, dtype=tables.dtype_from_descr(ProteinCacheInfo))
            h5.create_table("/", "ortholog_counts", ProteinCacheInfo, obj=cnts)
            logger.info("finished writing output file")


def combine_results(job_results, out):
    with tables.open_file(out, "w", filters=tables.Filters(complib="blosc2", complevel=5, fletcher32=True)) as fout:
        json_buffer = fout.create_earray(
            "/RootHOG", "JsonBuffer", tables.StringAtom(1), (0,), createparents=True, expectedrows=1e9
        )
        cnts, offsets, cur_off = [], [], 0
        for fn in job_results:
            with tables.open_file(fn, "r") as fin:
                json_buffer.append(fin.get_node("/family_json/buffer").read())
                off = fin.get_node("/family_json/offset").read()
                off["offset"] += cur_off
                cur_off += len(fin.get_node("/family_json/offset"))
                offsets.append(off)
                cnts.append(fin.get_node("/ortholog_counts").read())

        off = numpy.concatenate(offsets)
        if off["length"].sum() != len(json_buffer):
            logger.error(
                "Cached json seems broken. inconsistent lengths: %d <--> %d", off["length"].sum(), len(json_buffer)
            )
            raise DBConsistencyError("Cached json seems broken")
        off.sort(order="Fam")
        rhog_meta = numpy.zeros(len(off), dtype=tables.dtype_from_descr(RootHOGMetaTable))
        rhog_meta["FamNr"] = off["Fam"]
        rhog_meta["FamDataJsonOffset"] = off["offset"]
        rhog_meta["FamDataJsonLength"] = off["length"]
        roothog_meta = fout.create_table("/RootHOG", "MetaData", RootHOGMetaTable, obj=rhog_meta)
        roothog_meta.colinstances["FamNr"].create_csindex()

        cnts = numpy.concatenate(cnts)
        cnts.sort(order="EntryNr")
        if cnts["EntryNr"][0] == 0:
            cnts = cnts[cnts["EntryNr"] > 0]
        if len(cnts) != cnts[-1]["EntryNr"]:
            logger.error("Cached orthologs seem not complete: %d <--> %d", len(cnts), cnts[-1]["EntryNr"])
            raise DBConsistencyError("Cached orthologs seem not complete")

        tab = fout.create_table("/Protein", "OrthologsCountCache", ProteinCacheInfo, createparents=True, obj=cnts)
        tab.colinstances["EntryNr"].create_csindex()
