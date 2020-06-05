import tables
import numpy
import numpy.lib.recfunctions
import multiprocessing as mp
import collections
import logging

from tqdm import tqdm

from .db import Database
from .tablefmt import ProteinCacheInfo
from os.path import commonprefix, split

logger = logging.getLogger(__name__)

Protein = collections.namedtuple("Protein", ("entry_nr", "hog_id", "group"))


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

        for fun, params in iter(self.in_queue.get, None):
            self.out_queue.put(getattr(self, fun)(*params))
        self.out_queue.put(None)
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


def build_cache(db_fpath, nr_procs=None):
    request_queue = mp.Queue()
    result_queue = mp.Queue()
    nr_procs = nr_procs if nr_procs else mp.cpu_count()

    db = Database(db_fpath)
    nr_fams = db.get_nr_toplevel_hogs()
    singletons = [
        (r["EntryNr"], r["OmaGroup"])
        for r in db.get_hdf5_handle()
        .get_node("/Protein/Entries")
        .where('OmaHOG == b""')
    ]
    db.close()

    workers = []
    for i in range(nr_procs):
        w = CacheBuilderWorker(db_fpath, request_queue, result_queue, daemon=True)
        w.start()
        workers.append(w)

    for fam in range(nr_fams):
        request_queue.put(("analyse_fam", (fam,)))
    for singleton in singletons:
        request_queue.put(("analyse_singleton", singleton))
    # Sentinel objects to allow clean shutdown: 1 per worker.
    for i in range(nr_procs):
        request_queue.put(None)

    finished = 0
    results = []
    for res in tqdm(iter(result_queue.get, None), total=nr_fams + len(singletons)):
        if res is None:
            finished += 1
            if finished == nr_procs:
                break
        else:
            results.append(res)

    for w in workers:
        w.join()

    ret = numpy.lib.recfunctions.stack_arrays(results, usemask=False)
    ret.sort(order="EntryNr")
    return ret


def compute_and_store_cached_data(db_fpath, cache_path, nr_procs=None, force=False):
    with tables.open_file(db_fpath, "a") as h5:
        try:
            n = h5.get_node(cache_path)
            if force:
                h5.remove_node(cache_path)
            else:
                return
        except tables.NoSuchNodeError:
            pass

    cache = build_cache(db_fpath, nr_procs=nr_procs)

    path, name = split(cache_path)
    with tables.open_file(db_fpath, "a") as h5:
        tab = h5.create_carray(
            path, name, ProteinCacheInfo, createparents=True, obj=cache
        )
        tab.colinstances["EntryNr"].create_csindex()
