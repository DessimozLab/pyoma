import collections
import itertools
import logging
import os
import pickle
import re
import signal
import time
from functools import partial
import random

from tqdm import tqdm

import numpy
import tables
from datasketch import MinHash, MinHashLSH, MinHashLSHForest, LeanMinHash
from mpire import WorkerPool

from .db import Database

logger = logging.getLogger(__name__)
MinHash256 = partial(MinHash, seed=1, num_perm=256)
LeanMinHash256 = partial(LeanMinHash, seed=1)


class HogHasher(object):
    def __init__(self, db: Database):
        self.db = db
        self.xrefs = db.id_mapper["XRef"]

    def analyze_fam(self, fam_nr):
        members = self.db.member_of_fam(fam_nr)
        minhashes = collections.defaultdict(MinHash256)
        member_ids = []
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
            member_ids.append((hog_id, prot_id))
            if time.time() - t0 > 10:
                _log = logger.info if (time.time() - t_start) > 120 else logger.debug
                _log(f"working since {time.time()-t_start}sec on fam {fam_nr}. Done {i} out of {len(members)} proteins")
                t0 = time.time()
        return minhashes, member_ids


class HogHashWriter:
    """Write-only: buffers and flushes minhash + member data into an HDF5 file."""

    def __init__(self, hash_file, mode="w"):
        if mode not in ("a", "w"):
            raise ValueError(f"HogHashWriter only supports 'w' and 'a' modes, got '{mode}'")
        if mode == "a" and os.path.exists(hash_file):
            h5 = tables.open_file(hash_file, mode="a")
        else:
            h5 = self._create_hash_file(hash_file)
        self.h5 = h5
        self.hashes: tables.EArray = self.h5.get_node("/hashes")
        self.hogids: tables.EArray = self.h5.get_node("/hogids")
        try:
            self.members_table: tables.Table = self.h5.get_node("/members")
        except tables.NoSuchNodeError:
            self.members_table = None
        self.hash_buffer = []
        self.hogid_buffer = []
        self.member_buffer = []

    @staticmethod
    def _create_hash_file(hash_file):
        filters = tables.Filters(complevel=5, complib="blosc", bitshuffle=True, fletcher32=False, shuffle=True)
        h5 = tables.open_file(hash_file, mode="w", filters=filters)
        h5.create_earray("/", "hashes", atom=tables.Int64Atom(), shape=(0, 256), expectedrows=1_000_000)
        h5.create_earray(
            "/",
            "hogids",
            atom=tables.StringAtom(itemsize=255),
            shape=(0,),
            expectedrows=1_000_000,
        )
        h5.create_table(
            "/",
            "members",
            description=numpy.dtype([("OmaHOG", "S255"), ("ProtId", "S255")]),
            expectedrows=10_000_000,
        )
        return h5

    def add_minhashes(self, it, members=None):
        for hogid, minhash in it:
            self.hash_buffer.append(minhash.digest())
            self.hogid_buffer.append(hogid)
        if members:
            self.member_buffer.extend(members)
        if len(self.hash_buffer) > 10000:
            self.flush()

    def flush(self):
        if self.hash_buffer:
            self.hashes.append(self.hash_buffer)
            self.hogids.append(self.hogid_buffer)
            self.hash_buffer.clear()
            self.hogid_buffer.clear()
        if self.member_buffer and self.members_table is not None:
            self.members_table.append(self.member_buffer)
            self.member_buffer.clear()
        self.h5.flush()

    def close(self):
        self.flush()
        self.h5.close()


class HogLSHIndex:
    """Read-only: loads a hash file, builds an in-memory LSH index, and answers similarity queries."""

    def __init__(self, hash_file, threshold=0.7):
        if not os.path.exists(hash_file):
            raise IOError(f'file "{hash_file}" does not exist')
        self.h5 = tables.open_file(hash_file, mode="r")
        self.hashes: tables.EArray = self.h5.get_node("/hashes")
        self.hogids: tables.EArray = self.h5.get_node("/hogids")
        try:
            self.members_table: tables.Table = self.h5.get_node("/members")
        except tables.NoSuchNodeError:
            self.members_table = None
        self.threshold = threshold
        self.lsh = None
        self.hogid2row = None

    def compute_lsh(self):
        lsh = MinHashLSH(threshold=self.threshold, num_perm=256)
        hog2row = {}
        for row, (hogid, hashvals) in enumerate(itertools.zip_longest(self.hogids, self.hashes)):
            hog2row[hogid] = row
            lsh.insert(hogid, LeanMinHash256(hashvalues=hashvals))
        self.lsh = lsh
        self.hogid2row = hog2row

    def load_members(self):
        """Load /members into sorted in-memory arrays for exact Jaccard queries.

        OmaHOG values are kept as a sorted S255 array; ProtIds are folded to
        uint64 by XOR-ing their 8-byte chunks (fast, vectorised, negligible
        collision rate for realistic protein-ID strings).
        """
        if self.members_table is None:
            raise ValueError("This index has no /members table; rebuild the hash file to include member data.")
        hogids = self.members_table.col("OmaHOG")  # S255 ndarray
        protids = self.members_table.col("ProtId")  # S255 ndarray
        order = numpy.argsort(hogids, kind="stable")
        self._member_hogids = hogids[order]
        # Fold each 255-byte ProtId to a single uint64 by XOR-ing 8-byte chunks.
        raw = protids[order].view(numpy.uint8).reshape(-1, 255)
        # Pad to next multiple of 8 (256 bytes) so the view works cleanly.
        padded = numpy.zeros((len(raw), 256), dtype=numpy.uint8)
        padded[:, :255] = raw
        chunks = padded.view(numpy.uint64)  # shape (N, 32)
        self._member_protids = numpy.bitwise_xor.reduce(chunks, axis=1)

    def members_for(self, hog_id):
        """Return the uint64-hashed ProtIds for all members of hog_id (incl. sub-HOGs).

        Uses a prefix range scan on the sorted _member_hogids array built by
        load_members().  For cross-database exact-Jaccard queries pass the
        result as query_protids to query().
        """
        if isinstance(hog_id, str):
            hog_id = hog_id.encode()
        stripped = hog_id.rstrip(b"\x00")
        upper = stripped[:-1] + bytes([stripped[-1] + 1])
        lo = numpy.searchsorted(self._member_hogids, stripped)
        hi = numpy.searchsorted(self._member_hogids, upper)
        return self._member_protids[lo:hi]

    def query(self, key, minhash):
        """Yield (query-key, target-key, jaccard) for each LSH candidate.

        Jaccard is the MinHash estimate (~±6% for 256 perms).  For exact values
        call members_for() on both sides and pass to the module-level jaccard().
        """
        for c in self.lsh.query(minhash):
            hashvals = self.hashes[self.hogid2row[c]]
            h = LeanMinHash256(hashvalues=hashvals)
            yield key, c, minhash.jaccard(h)

    def close(self):
        self.h5.close()


class HogLSHForestIndex(HogLSHIndex):
    """Like HogLSHIndex but backed by MinHashLSHForest (no threshold, top-k queries)."""

    def compute_lsh(self):
        forest = MinHashLSHForest(num_perm=256)
        hog2row = {}
        for row, (hogid, hashvals) in enumerate(itertools.zip_longest(self.hogids, self.hashes)):
            hog2row[hogid] = row
            forest.add(hogid, LeanMinHash256(hashvalues=hashvals))
        forest.index()
        self.lsh = forest
        self.hogid2row = hog2row

    def query(self, key, minhash):
        for c in self.lsh.query(minhash, k=10):
            hashvals = self.hashes[self.hogid2row[c]]
            h = LeanMinHash256(hashvalues=hashvals)
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

    def _sort_by_size_descending(db_path, fams):
        with Database(db_path) as db:
            try:
                family_idx = db.get_hdf5_handle().get_node("/Protein/FamIndex/Lookup")
                sizes = family_idx.col("count")
                return sorted(fams, key=lambda f: sizes[f], reverse=True)
            except tables.NoSuchNodeError:
                logger.warning("FamIndex not available; processing families in random order")
                fams = list(fams)
                random.shuffle(fams)
        return list(fams)

    if lsh_path is None or not os.path.exists(lsh_path):
        fams_to_process = range(1, _get_nr_families(db_path) + 1)
    else:
        fams_to_process = _load_unprocessed_fams(db_path, lsh_path)
    return _sort_by_size_descending(db_path, fams_to_process)


class FamilyTimeoutError(Exception):
    pass


def hasher_worker_init(worker_state, db_path, log_conf, per_fam_timeout=1800):
    logging.basicConfig(level=log_conf["level"] + 5, format=log_conf["format"], datefmt=log_conf["datefmt"])
    worker_state["db_path"] = db_path
    worker_state["db"] = Database(db_path)
    worker_state["hasher"] = HogHasher(worker_state["db"])
    worker_state["per_fam_timeout"] = per_fam_timeout
    logger.info("initializing hasher worker %s", worker_state)


def _reinit_worker_db(worker_state):
    """Close and reopen the database after a timeout interruption.

    SIGALRM can leave the HDF5 file handle in an inconsistent state, so we
    must replace it before the worker picks up the next family.
    """
    try:
        worker_state["db"].close()
    except Exception:
        pass
    worker_state["db"] = Database(worker_state["db_path"])
    worker_state["hasher"] = HogHasher(worker_state["db"])
    logger.info("Worker database connection reinitialized after timeout")


def hash_worker_exit(worker_state):
    logger.info("exiting hash worker, closing database handle")
    del worker_state["hasher"]
    worker_state["db"].close()
    del worker_state["db"]


def hash_worker_fn(worker_state, fam):
    logger.debug("computing hashes for family %s", fam)
    t0 = time.time()
    timeout = worker_state.get("per_fam_timeout", 1200)

    def _timeout_handler(signum, frame):
        raise FamilyTimeoutError(fam)

    old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(timeout)
    try:
        hashes, members = worker_state["hasher"].analyze_fam(fam)
        signal.alarm(0)
        logger.debug("... done with fam %s. Took %f sec", fam, time.time() - t0)
        return hashes, members
    except FamilyTimeoutError:
        logger.warning("Skipping family %s: still running after %ds (likely PyTables index busy-loop)", fam, timeout)
        _reinit_worker_db(worker_state)
        return {}, []
    except Exception:
        signal.alarm(0)
        logger.exception("Error processing family %s", fam)
        raise
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


def get_logging_config():
    root = logging.getLogger()
    handler = next(h for h in root.handlers if h.formatter is not None)
    log_config = {
        "level": root.level,
        "format": handler.formatter._fmt,
        "datefmt": handler.formatter.datefmt,
    }
    return log_config


def compute_minhashes_for_db(db_path, output_path, nr_procs=None, per_fam_timeout=600):
    """Compute minhashes for all HOG families in db_path.

    per_fam_timeout: seconds before a stuck worker is interrupted and the family skipped.
                     Addresses the PyTables busy-loop bug on databases without FamIndex.
                     SIGALRM interrupts Python-level loops; C-level loops inside HDF5/numexpr
                     are not interruptible this way and would require a SIGKILL watchdog.
    """
    fams_to_process = generator_of_unprocessed_fams(db_path, output_path)
    collector = HogHashWriter(output_path, mode="a")
    try:
        with WorkerPool(n_jobs=nr_procs, use_worker_state=True, keep_alive=True, start_method="spawn") as pool:
            worker_init_ = partial(
                hasher_worker_init,
                db_path=db_path,
                log_conf=get_logging_config(),
                per_fam_timeout=per_fam_timeout,
            )
            results = pool.imap_unordered(
                hash_worker_fn,
                fams_to_process,
                worker_init=worker_init_,
                worker_exit=hash_worker_exit,
                worker_lifespan=10000,
                chunk_size=10,
                progress_bar=True,
            )
            for hashes, members in results:
                if hashes:
                    collector.add_minhashes(hashes.items(), members)
            collector.flush()
    finally:
        collector.close()


def jaccard(query_proteins: numpy.ndarray, cand_proteins: numpy.ndarray) -> float:
    qp = numpy.unique(query_proteins)
    cp = numpy.unique(cand_proteins)
    n_common = len(numpy.intersect1d(qp, cp, assume_unique=True))
    n_union = len(qp) + len(cp) - n_common
    return n_common / n_union if n_union > 0 else 0.0


def compare_versions(output_file, target_path, *old_path, index_cls=HogLSHIndex):
    lsh = index_cls(target_path)
    lsh.compute_lsh()
    use_exact = lsh.members_table is not None
    if use_exact:
        logger.info("Loading target member sets for exact Jaccard computation")
        lsh.load_members()
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
        for old_fn in old_path:
            old = index_cls(old_fn)
            if use_exact and old.members_table is not None:
                old.load_members()
            else:
                use_exact = False
            buffer = []
            c = collections.Counter()

            logger.info(
                f"{old_fn} --> |hogids|={len(old.hogids)}, |hashes|={len(old.hashes)}, use_exact={bool(use_exact)}"
            )
            for old_id, old_hashvals in tqdm(itertools.zip_longest(old.hogids, old.hashes), total=len(old.hogids)):
                minhash = LeanMinHash256(hashvalues=old_hashvals)
                raw_candidates = sorted(lsh.query(old_id, minhash), key=lambda x: -x[2])
                if raw_candidates:
                    raw_candidates = [c for c in raw_candidates if c[2] > 0.25 * raw_candidates[0][2]]
                if use_exact:
                    old_members = old.members_for(old_id)
                    candidates = sorted(
                        [(key, c, jaccard(old_members, lsh.members_for(c))) for key, c, _ in raw_candidates],
                        key=lambda x: -x[2],
                    )
                else:
                    candidates = raw_candidates
                c["query"] += 1
                if len(candidates) > 0 and candidates[0][1] != raw_candidates[0][1]:
                    logger.debug(
                        f"exact jaccard changed top hit: old_id={old_id.decode()}, min-hash-jaccard: {raw_candidates[0]}, exact-jaccard: {candidates[0]}"
                    )
                    c["invert"] += 1

                if len(candidates) > 10 and candidates[6][2] > 0.9:
                    # this is a dubious node, store it for now, maybe try to recover later.
                    # it seems that this happens if the CanonicalId is truncated and hence maps to several ids.
                    dubious.append([old_id])
                    c["dubious"] += 1
                    continue
                if len(candidates) > 0 and candidates[0][2] > 0.3:
                    buffer.append((old_id, candidates[0][1], candidates[0][2]))
                    c["mapped"] += 1
                    if candidates[0][2] < 1:
                        top_jaccard = candidates[0][2]
                        minor_matches = [(old_id, c[1], c[2]) for c in candidates[1:] if c[2] > 0.5 * top_jaccard]
                        if minor_matches:
                            buffer.extend(minor_matches)
                            c["has_minor"] += 1
                            c["minor_matches"] += len(minor_matches)
                else:
                    c["not_mapped"] += 1
                if len(buffer) > 10_000:
                    tab.append(buffer)
                    buffer.clear()
            if len(buffer) > 0:
                tab.append(buffer)
                buffer.clear()
            tab.flush()
            dubious.flush()
            old.close()

            msg = [f"{key:20s}: {val:7d}" for key, val in c.items()]
            logger.info(f"Mapping statistics:\n" + ("*" * 20) + "\n" + "\n".join(msg))
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
