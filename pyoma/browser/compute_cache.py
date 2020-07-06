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

from tqdm import tqdm

from .db import Database
from .tablefmt import ProteinCacheInfo
from os.path import commonprefix, split

logger = logging.getLogger(__name__)

Protein = collections.namedtuple("Protein", ("entry_nr", "hog_id", "group"))


def signal_handler(signum, frame):
    logger.info("received signal " + str(signum))
    raise KeyboardInterrupt(signum)


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

        try:
            for job in iter(self.in_queue.get, None):
                fun, params = job
                res = getattr(self, fun)(*params)
                self.out_queue.put((job, res))
            self.out_queue.put(None)
        except KeyboardInterrupt:
            logger.info("received interrupt. Terminating")
        self.db.close()

    def load_fam_members(self, fam):
        members = []
        hog_range = ["HOG:{:07d}".format(x).encode("utf8") for x in (fam, fam + 1)]
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
        grp_members = {
            grp: set(self.load_grp_members(grp))
            for grp in set(z.group for z in fam_members if z.group > 0)
        }
        counts = numpy.zeros(
            len(fam_members), dtype=tables.dtype_from_descr(ProteinCacheInfo)
        )
        for i, p1 in enumerate(fam_members):
            vps = set(self.load_vps(p1.entry_nr))
            ind_orth = set(p2.entry_nr for p2 in fam_members if are_orthologous(p1, p2))
            grp = grp_members.get(p1.group, set([])) - set([p1.entry_nr])
            logger.debug(
                "entry {}: vps: {} ipw: {} grp: {} any: {}".format(
                    p1.entry_nr, vps, ind_orth, grp, vps | ind_orth | grp
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


def build_cache(db_fpath, nr_procs=None, from_cache=None):
    request_queue = mp.Queue()
    result_queue = mp.Queue()
    nr_procs = nr_procs if nr_procs else mp.cpu_count()

    if from_cache is not None and os.path.isfile(from_cache):
        with tables.open_file(from_cache, "r") as cache:
            results = [cache.root.results.read()]
            jobs = list(cache.root.pending_jobs.read()[0])
        shutil.copy2(from_cache, from_cache + ".restarted")
        logger.info(
            "loaded results for {} proteins and {} remaining jobs".format(
                len(results[0]), len(jobs)
            )
        )
    else:
        db = Database(db_fpath)
        nr_fams = db.get_nr_toplevel_hogs()
        singletons = [
            (int(r["EntryNr"]), int(r["OmaGroup"]))
            for r in db.get_hdf5_handle()
            .get_node("/Protein/Entries")
            .where('OmaHOG == b""')
        ]
        nr_entries = len(db.get_hdf5_handle().get_node("/Protein/Entries"))
        db.close()

        logger.info(
            "found {} hog and {} singleton jobs to be computed".format(
                nr_fams, len(singletons)
            )
        )
        jobs = [("analyse_fam", (fam + 1,)) for fam in range(nr_fams)]
        jobs.extend([("analyse_singleton", singleton) for singleton in singletons])
        logger.info(
            "nr of jobs: {} (expected {})".format(len(jobs), nr_fams + len(singletons))
        )
        results = []

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
            if reply is None:
                finished += 1
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
                results.append(res)
                pending_jobs.remove(job)
                if time.time() - last_cache_timestamp > 60 and from_cache is not None:
                    sofar_results = numpy.lib.recfunctions.stack_arrays(
                        results, usemask=False
                    )
                    write_cache(from_cache, sofar_results, pending_jobs)
                    last_cache_timestamp = time.time()
        for w in workers:
            w.join()
    except KeyboardInterrupt as e:
        logger.info("recived interrupt. writeing out temp results")
        sofar_results = numpy.lib.recfunctions.stack_arrays(results, usemask=False)
        if from_cache is None:
            from_cache = "partial_compute_cache.h5"
        write_cache(from_cache, sofar_results, pending_jobs)
        sys.exit(99)

    ret = numpy.lib.recfunctions.stack_arrays(results, usemask=False)
    ret.sort(order="EntryNr")
    print(ret)
    assert check_all_there(nr_entries, ret)
    return ret


def write_cache(fn, sofar_results, pending_jobs):
    if os.path.exists(fn):
        os.replace(fn, fn + ".0")
    with tables.open_file(fn, "w") as h5:
        h5.create_table("/", "results", ProteinCacheInfo, obj=sofar_results)
        a = h5.create_vlarray("/", "pending_jobs", tables.ObjectAtom())
        a.append(pending_jobs)
        h5.flush()
    logger.info("written to cache {}".format(fn))


def check_all_there(nr_prots, cache):
    if len(cache) == nr_prots:
        return True
    missings = set(range(1, nr_prots + 1)) - set(cache["EntryNr"])
    logger.error("Missing cache value for {}".format(missings))
    return False


def compute_and_store_cached_data(
    db_fpath, cache_path, nr_procs=None, force=False, tmp_cache=None
):
    with tables.open_file(db_fpath, "a") as h5:
        try:
            n = h5.get_node(cache_path)
            if force:
                h5.remove_node(cache_path)
            else:
                return
        except tables.NoSuchNodeError:
            pass

    signal.signal(signal.SIGUSR2, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    cache = build_cache(db_fpath, nr_procs=nr_procs, from_cache=tmp_cache)

    path, name = split(cache_path)
    with tables.open_file(db_fpath, "a") as h5:
        tab = h5.create_table(
            path, name, ProteinCacheInfo, createparents=True, obj=cache
        )
        tab.colinstances["EntryNr"].create_csindex()
