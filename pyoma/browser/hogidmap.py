import collections
import inspect
import itertools
import logging
import os
import pickle
import re
import time
from functools import partial
from tqdm import tqdm

import numpy
import tables
from datasketch import MinHash, MinHashLSH, LeanMinHash
from mpire import WorkerPool

from .db import Database

logger = logging.getLogger(__name__)
MinHash256 = partial(MinHash, seed=1, num_perm=256)
LeanMinHash256 = partial(LeanMinHash, seed=1)

# datasketch >= 2.0 introduced a `scheme` parameter/attribute on MinHash and
# LeanMinHash (e.g. "affine32", or "legacy" for hash values produced by
# datasketch < 2.0, which predates the concept entirely). Reconstructing a
# MinHash/LeanMinHash from raw stored hashvalues now requires that scheme to
# be passed explicitly; datasketch < 2.0 has no such parameter and raises
# TypeError if given one.
_MINHASH_HAS_SCHEME = "scheme" in inspect.signature(MinHash.__init__).parameters


def _scheme_kwargs(scheme):
    """kwargs needed to rebuild a MinHash/LeanMinHash from stored hashvalues
    with the permutation `scheme` they were created with, compatible with
    both datasketch < 2.0 and >= 2.0."""
    if not _MINHASH_HAS_SCHEME:
        return {}
    return {"scheme": scheme or "legacy"}


class HogHasher(object):
    def __init__(self, db: Database):
        self.db = db
        self.xrefs = db.id_mapper["XRef"]

    def analyze_fam(self, fam_nr):
        members = self.db.member_of_fam(fam_nr)
        minhashes = collections.defaultdict(MinHash256)
        t_start = t0 = time.time()
        for i, e in enumerate(members):
            hog_id = e["OmaHOG"]
            prot_id = e["CanonicalId"]
            if len(prot_id) == e.dtype["CanonicalId"].itemsize:
                # probably overflow, get full id from xref
                for ref in self.xrefs.iter_xrefs_for_entry_nr(e["EntryNr"]):
                    if ref["xref"].encode("utf-8").startswith(prot_id):
                        prot_id = ref["xref"].encode("utf-8")
                        break
            for p in re.finditer(rb"\.", hog_id):
                minhashes[hog_id[: p.start()]].update(prot_id)
            minhashes[hog_id].update(prot_id)
            if time.time() - t0 > 10:
                _log = logger.info if (time.time() - t_start) > 120 else logger.debug
                _log(f"working since {time.time()-t_start}sec on fam {fam_nr}. Done {i} out of {len(members)} proteins")
                t0 = time.time()
        return minhashes


class LSHBuilder(object):
    def __init__(self, hash_file, mode="r", threshold=0.7):
        if mode not in ("r", "a", "w"):
            raise ValueError("invalid mode string ``%s``. Allowed modes are: " "'r', 'a' and 'w'" % mode)
        self._min_hash_scheme = None
        if mode == "r":
            if not os.path.exists(hash_file):
                raise IOError('file "{}" does not exist'.format(hash_file))
            self.h5, self.lsh, self.hogid2row = self._load_hash_file(hash_file, mode)
        elif mode == "a" and os.path.exists(hash_file):
            self.h5, self.lsh, self.hogid2row = self._load_hash_file(hash_file, mode)
        elif mode == "w" or mode == "a" and not os.path.exists(hash_file):
            self.lsh = MinHashLSH(threshold=threshold, num_perm=256)
            self.hogid2row = {}
            self.h5 = self.init_hash_table_file(hash_file)
        self.hashes: tables.EArray = self.h5.get_node("/hashes")
        self.hogids: tables.EArray = self.h5.get_node("/hogids")
        self.threshold = threshold
        self._readonly = mode == "r"

    def _open_hdf5(self, filename, mode="w"):
        filters = None
        if mode == "w":
            filters = tables.Filters(complevel=5, complib="blosc", bitshuffle=True, fletcher32=False, shuffle=True)
        return tables.open_file(filename, mode=mode, filters=filters)

    def init_hash_table_file(self, hash_file):
        h5 = self._open_hdf5(hash_file, mode="w")
        h5.create_earray("/", "hashes", atom=tables.Int64Atom(), shape=(0, 256), expectedrows=1_000_000)
        h5.create_earray(
            "/",
            "hogids",
            atom=tables.StringAtom(itemsize=255),
            shape=(0,),
            expectedrows=1_000_000,
        )
        # Use VLArray of UInt8Atom to store serialized LSH object
        h5.create_vlarray("/", "lsh_obj", atom=tables.UInt8Atom())
        return h5

    def _load_hash_file(self, hash_file, mode="r"):
        h5 = self._open_hdf5(hash_file, mode=mode)
        lsh_obj_arr: tables.VLArray = h5.get_node("/lsh_obj")
        lsh_bytes = bytes(lsh_obj_arr[0])
        hogid2row_bytes = bytes(lsh_obj_arr[1])
        lsh = pickle.loads(lsh_bytes)
        hogid2row = pickle.loads(hogid2row_bytes)
        # Files written before this attribute existed predate datasketch's
        # scheme concept entirely, i.e. they were produced with the "legacy"
        # permutation scheme.
        self._min_hash_scheme = getattr(h5.get_node("/hashes")._v_attrs, "minhash_scheme", "legacy")
        return h5, lsh, hogid2row

    def close(self):
        if self.h5.mode != "r":
            self.hashes._v_attrs.minhash_scheme = self._min_hash_scheme or "legacy"
            lsh_obj_arr: tables.VLArray = self.h5.get_node("/lsh_obj")
            # Remove old content, if any (for repeated closes)
            lsh_obj_arr.truncate(0)
            # Serialize LSH and hogid2row
            lsh_bytes = pickle.dumps(self.lsh)
            hogid2row_bytes = pickle.dumps(self.hogid2row)
            lsh_obj_arr.append(numpy.frombuffer(lsh_bytes, dtype=numpy.uint8))
            lsh_obj_arr.append(numpy.frombuffer(hogid2row_bytes, dtype=numpy.uint8))
            lsh_obj_arr.flush()
            self.hashes.flush()
            self.hogids.flush()
        self.h5.close()

    def add_minhashes(self, it):
        hash_buffer = []
        hogid_buffer = []
        for hogid, minhash in it:
            scheme = getattr(minhash, "scheme", "legacy")
            if self._min_hash_scheme is None:
                self._min_hash_scheme = scheme
            elif self._min_hash_scheme != scheme:
                raise ValueError(
                    "MinHash scheme %r does not match scheme %r of previously added hashes in this "
                    "file; cannot mix minhashes computed with different datasketch schemes"
                    % (scheme, self._min_hash_scheme)
                )
            hash_buffer.append(minhash.digest())
            hogid_buffer.append(hogid)
        self.hashes.append(hash_buffer)
        self.hogids.append(hogid_buffer)
        self.hashes.flush()
        self.hogids.flush()

    def compute_lsh(self):
        lsh = MinHashLSH(threshold=self.threshold, num_perm=256)
        hog2row = {}
        for row, (hogid, hashvals) in enumerate(itertools.zip_longest(self.hogids, self.hashes)):
            hog2row[hogid] = row
            lsh.insert(hogid, LeanMinHash256(hashvalues=hashvals, **_scheme_kwargs(self._min_hash_scheme)))
        self.lsh = lsh
        self.hogid2row = hog2row

    def query(self, key, minhash):
        """query with a minhash, returns a list of tuples of
        (query-key, target-key, jaccard)"""
        candidates = self.lsh.query(minhash)
        for c in candidates:
            hashvals = self.hashes[self.hogid2row[c]]
            h = MinHash256(hashvalues=hashvals, **_scheme_kwargs(self._min_hash_scheme))
            yield key, c, minhash.jaccard(h)


def generator_of_unprocessed_fams(db_path, lsh_path=None):
    def _get_nr_families(db_path):
        with Database(db_path) as db:
            nr_hogs = db.get_nr_toplevel_hogs()
        logger.info("Found %d families to process", nr_hogs)
        return nr_hogs

    def _load_unprocessed_fams(db_path, lsh_path):
        with Database(db_path) as db:
            with tables.open_file(lsh_path, "r") as lsh_h5:
                processed_fams = set(db.parse_hog_id(x) for x in lsh_h5.get_node("/hogids"))
            remaining = set(range(1, db.get_nr_toplevel_hogs() + 1)) - processed_fams
        logger.info("Found %d unprocessed families", len(remaining))
        return remaining

    if lsh_path is None or not os.path.exists(lsh_path):
        fams_to_process = range(1, _get_nr_families(db_path) + 1)
    else:
        fams_to_process = _load_unprocessed_fams(db_path, lsh_path)
    return fams_to_process


def hasher_worker_init(worker_state, db_path, log_conf):
    logging.basicConfig(level=log_conf["level"] + 5, format=log_conf["format"], datefmt=log_conf["datefmt"])
    worker_state["db"] = Database(db_path)
    worker_state["hasher"] = HogHasher(worker_state["db"])
    logger.info("initializing hasher worker %s", worker_state)


def hash_worker_exit(worker_state):
    logger.info("exiting hash worker, closing database handle")
    del worker_state["hasher"]
    worker_state["db"].close()
    del worker_state["db"]


def hash_worker_fn(worker_state, fam):
    logger.debug("computing hashes for family %s", fam)
    t0 = time.time()
    try:
        hashes = worker_state["hasher"].analyze_fam(fam)
        logger.debug("... done with fam %s. Took %f sec", fam, time.time() - t0)
        return hashes
    except Exception:
        logger.exception("Error processing family %s", fam)
        raise


def get_logging_config():
    root = logging.getLogger()
    handler = next(h for h in root.handlers if h.formatter is not None)
    log_config = {
        "level": root.level,
        "format": handler.formatter._fmt,
        "datefmt": handler.formatter.datefmt,
    }
    return log_config


def compute_minhashes_for_db(db_path, output_path, nr_procs=None):
    fams_to_process = generator_of_unprocessed_fams(db_path, output_path)
    collector = LSHBuilder(output_path, mode="a")
    with WorkerPool(n_jobs=nr_procs, use_worker_state=True, keep_alive=True, start_method="spawn") as pool:
        worker_init_ = partial(hasher_worker_init, db_path=db_path, log_conf=get_logging_config())
        results = pool.imap_unordered(
            hash_worker_fn,
            fams_to_process,
            worker_init=worker_init_,
            worker_exit=hash_worker_exit,
            worker_lifespan=10000,
            chunk_size=100,
            progress_bar=True,
        )
        for hashes in results:
            collector.add_minhashes(hashes.items())
    collector.close()

    # while True:
    #     pipeline = Pipeline()
    #     pipeline.add_stage(Stage(FamGenerator, nr_procs=1, fam_generator=fams_to_process))
    #     pipeline.add_stage(Stage(HashWorker, nr_procs=nr_procs, db_path=db_path))
    #     pipeline.add_stage(Stage(Collector, nr_procs=1, output_path=output_path))
    #     print("setup pipeline, about to start it.")
    #     pipeline.run()
    #     print("finished with computing the MinHashLSH for {}".format(db_path))
    #
    #     fams_to_process = set(generator_of_unprocessed_fams(db_path, output_path))
    #     if len(fams_to_process) == 0:
    #         break
    # print("for sure all families processes. we're done!")


def compare_versions(output_file, target_path, *old_path):
    lsh = LSHBuilder(target_path, mode="r")
    lsh.compute_lsh()
    with tables.open_file(
        output_file,
        "w",
        filters=tables.Filters(complevel=7, complib="blosc2", shuffle=True),
    ) as h5_map:
        tab = h5_map.create_table(
            "/",
            "hogmap",
            description=numpy.dtype([("Old", "S255"), ("New", "S255"), ("Jaccard", "f4")]),
        )
        dubious = h5_map.create_earray(
            "/",
            "dubious",
            atom=tables.StringAtom(itemsize=255),
            shape=(0,),
            expectedrows=1e5,
        )
        for old in old_path:
            old = LSHBuilder(old, mode="r")
            for old_id, old_hashvals in tqdm(itertools.zip_longest(old.hogids, old.hashes), total=len(old.hogids)):
                minhash = LeanMinHash256(hashvalues=old_hashvals, **_scheme_kwargs(old._min_hash_scheme))
                candidates = sorted(lsh.query(old_id, minhash), key=lambda x: -x[2])
                logger.debug("old_id: %s: candidates: %s", old_id, candidates)
                if len(candidates) > 10 and candidates[6][2] > 0.9:
                    # this is a dubious node, store it for now, maybe try to recover later.
                    # it seems that this happens if the CanonicalId is truncated and hence maps to several ids.
                    dubious.append([old_id])
                    continue
                if len(candidates) > 0 and candidates[0][2] > 0.6:
                    tab.append([(old_id, candidates[0][1], candidates[0][2])])
                    if candidates[0][2] < 1:
                        other_cands = [(old_id, c[1], c[2]) for c in candidates[1:]]
                        if len(other_cands) > 0:
                            tab.append(other_cands)
            tab.flush()
            dubious.flush()
            old.close()
        # build index of Old ids
        tab.colinstances["Old"].create_csindex()


def build_lookup(target_db, old_dbs, nr_procs=None):
    def lsh_fn_from_db_path(dbpath, modif=None):
        dbname = os.path.splitext(os.path.basename(dbpath))[0]
        modif_str = "-{}".format(modif) if modif is not None else ""
        lsh_fn = "{}{}.hog-lsh.h5".format(dbname, modif_str)
        return lsh_fn

    target_lsh_fn = lsh_fn_from_db_path(target_db)
    compute_minhashes_for_db(target_db, target_lsh_fn, nr_procs=nr_procs)

    old_lsh_paths = []
    for k, dbpath in enumerate(old_dbs):
        cached_lsh_fn = dbpath.replace(".h5", ".hog-lsh.h5")
        if not os.path.exists(cached_lsh_fn):
            cached_lsh_fn = lsh_fn_from_db_path(dbpath, modif=k)
            print("computing lsh for {} - storing in {}".format(dbpath, cached_lsh_fn))
            compute_minhashes_for_db(dbpath, cached_lsh_fn, nr_procs=nr_procs)
        old_lsh_paths.append(cached_lsh_fn)

    hogmap_name = os.path.splitext(os.path.basename(target_db))[0] + ".hogmap.h5"
    compare_versions(hogmap_name, target_lsh_fn, old_lsh_paths)
