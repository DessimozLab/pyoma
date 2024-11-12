import collections
import logging
import itertools
import os
import pickle
import re
import json
from time import time

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
        for row in h5.get_node("/Protein/Entries"):
            if len(row["OmaHOG"]) > 0:
                fam_sizes.update((int(re_fam.match(row["OmaHOG"]).group(1)),))
            else:
                singletons.append(
                    (
                        int(row["EntryNr"]),
                        int(row["OmaGroup"]),
                    )
                )
    with open(out_prefix + "_singleton.pkl", "wb") as fh:
        pickle.dump(["process_singletons", singletons], fh)

    def yield_buckets(counts, nr_elem=10000):
        cur_bucket, cur_size = [], 0
        target = nr_elem**2
        for fam, size in sorted(counts.items(), key=lambda x: -x[1]):
            if size**2 > target:
                for rng in range(0, size, nr_elem):
                    yield [(fam, (rng, rng + nr_elem))]
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


def process_job_file(job_file: os.PathLike, db_fpath: os.PathLike, out: os.PathLike):
    with open(job_file, "rb") as fh:
        jobdata = pickle.load(fh)
    job, payload = jobdata
    with CacheBuilder(db_fpath, out) as builder:
        func = getattr(builder, job)
        if job == "process_singletons":
            func(payload)
        else:
            for args in payload:
                func(*args)


class CacheBuilder:
    def __init__(self, db_fpath, out_path):
        self.db_fpath = db_fpath
        self.out_path = out_path
        self.db = None
        self.h5 = None
        self.cnts = []
        self.json_buffer = []
        self.json_offsets = []
        self._buffer_offset = 0
        self._offset_dtype = [("Fam", "i4"), ("offset", "i8"), ("length", "i4")]

    def __enter__(self):
        self.db = Database(self.db_fpath)
        self.h5 = self.db.get_hdf5_handle()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.db.close()
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

        counts = numpy.zeros(nr_memb, dtype=tables.dtype_from_descr(ProteinCacheInfo))
        for i, p1 in tqdm(
            enumerate(fam_iter),
            disable=len(fam_members) < 500,
            desc=f"Processing family {fam}",
            total=nr_memb,
        ):
            vps = set(self.load_vps(p1.entry_nr))
            ind_orth = set(p2.entry_nr for p2 in fam_members if are_orthologous(p1, p2))
            grp = grp_members.get(p1.group, set([])) - set([p1.entry_nr])
            counts[i]["EntryNr"] = p1.entry_nr
            counts[i]["NrPairwiseOrthologs"] = len(vps)
            counts[i]["NrHogInducedPWOrthologs"] = len(ind_orth)
            counts[i]["NrHogInducedPWParalogs"] = len(fam_members) - len(ind_orth) - 1
            counts[i]["NrOMAGroupOrthologs"] = len(grp)
            counts[i]["NrAnyOrthologs"] = len(vps | ind_orth | grp)
        return counts

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
        roothog_meta.colinstances["Fam"].create_csindex()

        cnts = numpy.concatenate(cnts)
        cnts.sort(order="EntryNr")
        if len(cnts) != cnts[-1]["EntryNr"]:
            logger.error("Cached orthologs seem not complete: %d <--> %d", len(cnts), cnts[-1]["EntryNr"])
            raise DBConsistencyError("Cached orthologs seem not complete")

        tab = fout.create_table("/Protein", "OrthologsCountCache", ProteinCacheInfo, createparents=True, obj=cnts)
        tab.colinstances["EntryNr"].create_csindex()
