import tables
import numpy
import numpy.lib.recfunctions
import multiprocessing as mp
import collections
import logging
import signal
import shutil
import time
import os
import sys
import functools
import json

from tqdm import tqdm
from queue import Empty
from .db import Database
from .models import ProteinEntry
from .tablefmt import ProteinCacheInfo
from os.path import commonprefix, split

logger = logging.getLogger(__name__)

Protein = collections.namedtuple("Protein", ("entry_nr", "hog_id", "group"))


def signal_handler(signum, frame):
    logger.info("received signal " + str(signum))
    raise KeyboardInterrupt(signum)


def length_limited_set_formatter(s):
    if len(s) > 10:
        x = s.pop()
        s.add(x)
        return "Set(size={}, ex: {})".format(len(s), x)
    else:
        return str(s)


def are_orthologous(a: Protein, b: Protein):
    if a.entry_nr == b.entry_nr:
        return False
    prefix = commonprefix((a.hog_id, b.hog_id))
    if "." in prefix and prefix[-1].isdigit():
        return False
    return True


class CacheBuilderWorker(mp.Process):
    def __init__(self, db_fpath, in_queue, out_queue, **kwargs):
        super(CacheBuilderWorker, self).__init__(**kwargs)
        self.db_fpath = db_fpath
        self.in_queue = in_queue
        self.out_queue = out_queue
        self.h5 = None
        self.db = None

    def run(self):
        self.db = Database(self.db_fpath)
        self.h5 = self.db.get_hdf5_handle()
        timelimit_get_from_in_queue = functools.partial(self.in_queue.get, timeout=120)

        try:
            try:
                for job in iter(timelimit_get_from_in_queue, None):
                    fun, params = job
                    res = getattr(self, fun)(*params)
                    logger.debug("result for job ({}) ready".format(job))
                    self.out_queue.put((job, res))
            except Empty:
                logger.warning(
                    "No item nor termination signal received in Queue. Giving up"
                )
                logger.exception("Work-queue is empty")
            self.out_queue.put("DONE")
        except KeyboardInterrupt:
            logger.info("received interrupt. Terminating")
        self.db.close()
        logger.info("terminating worker process {}".format(self.name))

    def load_fam_members(self, fam):
        members = []
        hog_range = [self.db.format_hogid(x).encode("utf-8") for x in (fam, fam + 1)]
        for row in self.h5.get_node("/Protein/Entries").where(
            "({!r} <= OmaHOG) & (OmaHOG < {!r})".format(*hog_range)
        ):
            members.append(
                Protein(row["EntryNr"], row["OmaHOG"].decode(), row["OmaGroup"])
            )
        return members

    def load_vps(self, entry_nr):
        return self.db.get_vpairs(entry_nr)["EntryNr2"]

    def load_grp_members(self, group):
        return [
            row["EntryNr"]
            for row in self.h5.get_node("/Protein/Entries").where(
                "OmaGroup == {:d}".format(group)
            )
        ]

    def analyse_fam(self, fam):
        logger.debug("analysing family {}".format(fam))
        fam_members = self.load_fam_members(fam)
        logger.debug("family {} with {} members".format(fam, len(fam_members)))
        grp_members = {
            grp: set(self.load_grp_members(grp))
            for grp in set(z.group for z in fam_members if z.group > 0)
        }
        counts = numpy.zeros(
            len(fam_members), dtype=tables.dtype_from_descr(ProteinCacheInfo)
        )
        for i, p1 in tqdm(
            enumerate(fam_members),
            disable=len(fam_members) < 500,
            desc="fam {}".format(fam),
        ):
            vps = set(self.load_vps(p1.entry_nr))
            ind_orth = set(p2.entry_nr for p2 in fam_members if are_orthologous(p1, p2))
            grp = grp_members.get(p1.group, set([])) - set([p1.entry_nr])
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug(
                    "entry {}: vps: {} ipw: {} grp: {} any: {}".format(
                        p1.entry_nr,
                        length_limited_set_formatter(vps),
                        length_limited_set_formatter(ind_orth),
                        length_limited_set_formatter(grp),
                        length_limited_set_formatter(vps | ind_orth | grp),
                    )
                )
            counts[i]["EntryNr"] = p1.entry_nr
            counts[i]["NrPairwiseOrthologs"] = len(vps)
            counts[i]["NrHogInducedPWOrthologs"] = len(ind_orth)
            counts[i]["NrHogInducedPWParalogs"] = len(fam_members) - len(ind_orth) - 1
            counts[i]["NrOMAGroupOrthologs"] = len(grp)
            counts[i]["NrAnyOrthologs"] = len(vps | ind_orth | grp)
        return counts

    def analyse_singleton(self, entry_nr, group_nr):
        logger.debug("analysing singleton {} (grp {})".format(entry_nr, group_nr))
        vps = set(self.load_vps(entry_nr))
        grp_members = set([])
        if group_nr > 0:
            grp_members = set(self.load_grp_members(group_nr))
        counts = numpy.array(
            [(entry_nr, len(vps), 0, 0, len(grp_members), len(vps | grp_members))],
            dtype=tables.dtype_from_descr(ProteinCacheInfo),
        )
        return counts

    def compute_familydata_json(self, fam):
        famhog_id = self.db.format_hogid(fam)
        logger.debug("analysing family {}".format(fam))
        fam_members = self.load_fam_members(fam)
        logger.debug("family {} with {} members".format(fam, len(fam_members)))
        if len(fam_members) > 50000:
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
                print("gene_similarity failed for {}: {}".format(fam, e))
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


class ConsistenceyError(Exception):
    pass


class ResultHandler:
    def __init__(self, cache_file):
        self.cache_file = cache_file
        self.ortholog_count_result = []
        self.family_json_offset = []
        self.in_memory_json_buffer = []
        self.buffer_offset = 0
        self.jobs = []
        self._last_cache_timestamp = time.time()
        self._offset_dtype = [("Fam", "i4"), ("offset", "i8"), ("length", "i4")]

    def load_cache(self):
        def resilient_load_data_node(h5, node):
            try:
                res = [h5.get_node(node).read()]
            except tables.NoSuchNodeError:
                res = []
            return res

        with tables.open_file(self.cache_file) as cache:
            self.ortholog_count_result = resilient_load_data_node(
                cache, "/ortholog_counts"
            )
            self.family_json_offset = resilient_load_data_node(
                cache, "/family_json/offset"
            )
            self.buffer_offset = len(cache.root.family_json.buffer)
            self.jobs = cache.root.pending_jobs.read(0)[0]

    def add_jobs(self, jobs):
        self.jobs.extend(jobs)

    def handle_result(self, job, result):
        if job[0] == "compute_familydata_json":
            self.store_familydata_json_result(job, result)
        elif job[0] in ("analyse_fam", "analyse_singleton"):
            self.store_ortholog_count_result(job, result)
        else:
            raise ValueError("Unexpected result type")
        self.jobs.remove(job)
        if (
            time.time() - self._last_cache_timestamp > 300
            or sum(len(z) for z in self.in_memory_json_buffer) > 100e6
        ):
            self.write_to_disk()

    def store_familydata_json_result(self, job, result):
        fam = job[1][0]
        encoded_json = result.encode("utf-8")
        json_as_np = numpy.ndarray(
            (len(encoded_json),), buffer=encoded_json, dtype=tables.StringAtom(1)
        )
        self.in_memory_json_buffer.append(json_as_np)
        self.family_json_offset.append(
            numpy.array(
                [(fam, self.buffer_offset, len(encoded_json))], dtype=self._offset_dtype
            )
        )
        self.buffer_offset += len(encoded_json)

    def store_ortholog_count_result(self, job, result):
        self.ortholog_count_result.append(result)

    def write_to_disk(self):
        logger.info("writing a milestone to disk...")
        transfer_data = False
        if os.path.exists(self.cache_file):
            os.replace(self.cache_file, self.cache_file + ".0")
            transfer_data = True
        with tables.open_file(self.cache_file, "w") as h5:
            buf = h5.create_earray(
                "/family_json",
                "buffer",
                tables.StringAtom(1),
                (0,),
                createparents=True,
                expectedrows=1e9,
            )
            if transfer_data:
                with tables.open_file(self.cache_file + ".0", "r") as prev:
                    buf.append(prev.root.family_json.buffer.read())
            for el in self.in_memory_json_buffer:
                buf.append(el)
            buf.flush()
            if len(self.family_json_offset) > 0:
                off = numpy.lib.recfunctions.stack_arrays(
                    self.family_json_offset, usemask=False
                )
                h5.create_table("/family_json", "offset", None, obj=off)

            if len(self.ortholog_count_result) > 0:
                cnts = numpy.lib.recfunctions.stack_arrays(
                    self.ortholog_count_result, usemask=False
                )
                h5.create_table("/", "ortholog_counts", ProteinCacheInfo, obj=cnts)

            a = h5.create_vlarray("/", "pending_jobs", tables.ObjectAtom())
            a.append(self.jobs)
            h5.flush()
            if len(buf) != self.buffer_offset:
                raise ConsistenceyError(
                    "buffer has unexpeced length: {}vs{}".format(
                        len(buf), self.buffer_offset
                    )
                )
            self.in_memory_json_buffer = []
            self._last_cache_timestamp = time.time()
            logger.info("finished writing milestone to {}".format(self.cache_file))


def build_cache(db_fpath, nr_procs=None, from_cache=None):
    request_queue = mp.Queue()
    result_queue = mp.Queue()
    nr_procs = nr_procs if nr_procs else mp.cpu_count()

    db = Database(db_fpath)
    nr_entries = len(db.get_hdf5_handle().get_node("/Protein/Entries"))
    result_handler = ResultHandler(from_cache)
    if from_cache is not None and os.path.isfile(from_cache):
        result_handler.load_cache()
        jobs = result_handler.jobs
        logger.debug(jobs)
    else:
        nr_fams = db.get_nr_toplevel_hogs()
        singletons = [
            (int(r["EntryNr"]), int(r["OmaGroup"]))
            for r in db.get_hdf5_handle()
            .get_node("/Protein/Entries")
            .where('OmaHOG == b""')
        ]

        logger.info(
            "found {} hog and {} singleton jobs to be computed".format(
                nr_fams, len(singletons)
            )
        )
        jobs = [("analyse_fam", (fam + 1,)) for fam in range(nr_fams)]
        jobs.extend([("compute_familydata_json", (fam + 1,)) for fam in range(nr_fams)])
        jobs.extend([("analyse_singleton", singleton) for singleton in singletons])
        logger.info(
            "nr of jobs: {} (expected {})".format(
                len(jobs), 2 * nr_fams + len(singletons)
            )
        )
        result_handler.add_jobs(jobs)
    db.close()

    workers = []
    for i in range(nr_procs):
        w = CacheBuilderWorker(db_fpath, request_queue, result_queue, daemon=True)
        w.start()
        workers.append(w)

    pending_jobs = set([])
    for job in jobs:
        request_queue.put(job)
        pending_jobs.add(job)
    # Sentinel objects to allow clean shutdown: 1 per worker.
    for i in range(nr_procs):
        request_queue.put(None)

    finished = 0
    try:
        logger.info("start to receive results")
        last_cache_timestamp = time.time()
        for reply in tqdm(iter(result_queue.get, None), total=len(jobs)):
            if reply == "DONE":
                finished += 1
                logger.info("{} workers finished".format(finished))
                if finished == nr_procs:
                    if len(pending_jobs) > 0:
                        logger.error(
                            "still {} pending jobs...: {}".format(
                                len(pending_jobs), pending_jobs
                            )
                        )
                    break
            else:
                job, res = reply
                logger.debug("received result for job {})".format(job))
                result_handler.handle_result(job, res)
        result_handler.write_to_disk()
        logger.info("exit receiver loop. joining workers...")
        for w in workers:
            w.join()
        logger.debug("all workers joined")
    except KeyboardInterrupt as e:
        logger.info("recived interrupt. writeing out temp results")
        result_handler.write_to_disk()
        sys.exit(99)

    ret = numpy.lib.recfunctions.stack_arrays(
        result_handler.ortholog_count_result, usemask=False
    )
    ret.sort(order="EntryNr")
    logger.info("sorted results: {}".format(ret))
    assert check_all_there(nr_entries, ret)
    return ret


def check_all_there(nr_prots, cache):
    if len(cache) == nr_prots:
        return True
    missings = set(range(1, nr_prots + 1)) - set(cache["EntryNr"])
    logger.error("Missing cache value for {}".format(missings))
    return False


def compute_and_store_cached_data(db_fpath, nr_procs=None, force=False, tmp_cache=None):
    ortholog_cnts_cache_path = "/Protein/OrthologsCountCache"
    if tmp_cache is None:
        tmp_cache = "/tmp/compute_cache.h5"
    with tables.open_file(db_fpath, "a") as h5:
        try:
            n = h5.get_node(ortholog_cnts_cache_path)
            if force:
                h5.remove_node(ortholog_cnts_cache_path)
            else:
                return
        except tables.NoSuchNodeError:
            pass

    signal.signal(signal.SIGUSR2, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    cache = build_cache(db_fpath, nr_procs=nr_procs, from_cache=tmp_cache)
    update_hdf5_from_cachefile(db_fpath, tmp_cache)


def update_hdf5_from_cachefile(db_fpath, tmp_cache):
    ortholog_cnts_cache_path = "/Protein/OrthologsCountCache"
    path, name = split(ortholog_cnts_cache_path)

    with tables.open_file(db_fpath, "a") as h5, tables.open_file(
        tmp_cache, "r"
    ) as cache_h5:
        cached_cnts = cache_h5.get_node("/ortholog_counts").read()
        cached_cnts.sort(order="EntryNr")

        tab = h5.create_table(
            path, name, ProteinCacheInfo, createparents=True, obj=cached_cnts
        )
        tab.colinstances["EntryNr"].create_csindex()

        json_off = cache_h5.get_node("/family_json/offset").read()
        json_off.sort(order="Fam")
        tab = h5.root.RootHOG.MetaData.read()
        for fam, off, length in json_off:
            if tab[fam - 1]["FamNr"] != fam:
                raise ConsistenceyError("table not properly ordered")
            tab[fam - 1]["FamDataJsonOffset"] = off
            tab[fam - 1]["FamDataJsonLength"] = length
        h5.root.RootHOG.MetaData.modify_column(
            column=tab["FamDataJsonOffset"], colname="FamDataJsonOffset"
        )
        h5.root.RootHOG.MetaData.modify_column(
            column=tab["FamDataJsonLength"], colname="FamDataJsonLength"
        )

        json_in_buf = cache_h5.root.family_json.buffer
        json_in_buf._f_copy(
            h5.root.RootHOG, "JsonBuffer", expectedrows=len(json_in_buf)
        )
