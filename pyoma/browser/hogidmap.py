import collections
import itertools
import logging
import multiprocessing
import os
import re
import time
from functools import partial
from tqdm import tqdm

import numpy
import tables
from datasketch import MinHash, MinHashLSH, LeanMinHash

from .db import Database
from .hogprofile.build import BaseProfileBuilderProcess, SourceProcess, Pipeline, Stage

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
        for e in members:
            hog_id = e["OmaHOG"]
            prot_id = e["CanonicalId"]
            if len(prot_id) == e.dtype["CanonicalId"].itemsize:
                # probably overflow, get full id from xref
                for ref in self.xrefs.iter_xrefs_for_entry_nr(e["EntryNr"]):
                    if ref["xref"].encode("utf-8").startswith(prot_id):
                        prot_id = ref["xref"].encode("utf-8")
                        break
            for p in re.finditer(br"\.", hog_id):
                minhashes[hog_id[: p.start()]].update(prot_id)
            minhashes[hog_id].update(prot_id)
        return minhashes


class LSHBuilder(object):
    def __init__(self, hash_file, mode="r", threshold=0.7):
        if mode not in ("r", "a", "w"):
            raise ValueError(
                "invalid mode string ``%s``. Allowed modes are: "
                "'r', 'a' and 'w'" % mode
            )
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
        self.hashes = self.h5.get_node("/hashes")
        self.hogids = self.h5.get_node("/hogids")
        self.threshold = threshold
        self._readonly = mode == "r"

    def _open_hdf5(self, filename, mode="w"):
        filters = None
        if mode == "w":
            filters = tables.Filters(
                complevel=1, complib="blosc", shuffle=True, fletcher32=True
            )
        return tables.open_file(filename, mode=mode, filters=filters)

    def init_hash_table_file(self, hash_file):
        h5 = self._open_hdf5(hash_file, mode="w")
        h5.create_earray(
            "/", "hashes", atom=tables.Int64Atom(), shape=(0, 256), expectedrows=1e6
        )
        h5.create_earray(
            "/",
            "hogids",
            atom=tables.StringAtom(itemsize=255),
            shape=(0,),
            expectedrows=1e6,
        )
        h5.create_vlarray("/", "lsh_obj", atom=tables.ObjectAtom())
        return h5

    def _load_hash_file(self, hash_file, mode="r"):
        h5 = self._open_hdf5(hash_file, mode=mode)
        lsh_obj_arr = h5.get_node("/lsh_obj")
        lsh = lsh_obj_arr[-2]
        hogid2row = lsh_obj_arr[-1]
        return h5, lsh, hogid2row

    def close(self):
        if self.h5.mode != "r":
            lsh_obj_arr = self.h5.get_node("/lsh_obj")
            lsh_obj_arr.append(self.lsh)
            lsh_obj_arr.append(self.hogid2row)
            self.hashes.flush()
            self.hogids.flush()
        self.h5.close()

    def add_minhashes(self, it):
        for hogid, minhash in it:
            self.hashes.append([minhash.digest()])
            self.hogids.append([hogid])
        self.hashes.flush()
        self.hogids.flush()

    def compute_lsh(self):
        lsh = MinHashLSH(threshold=self.threshold, num_perm=256)
        hog2row = {}
        for row, (hogid, hashvals) in enumerate(
            itertools.zip_longest(self.hogids, self.hashes)
        ):
            hog2row[hogid] = row
            lsh.insert(hogid, LeanMinHash256(hashvalues=hashvals))
        self.lsh = lsh
        self.hogid2row = hog2row

    def query(self, key, minhash):
        """query with a minhash, returns a list of tuples of
        (query-key, target-key, jaccard)"""
        candidates = self.lsh.query(minhash)
        for c in candidates:
            hashvals = self.hashes[self.hogid2row[c]]
            h = MinHash256(hashvalues=hashvals)
            yield key, c, minhash.jaccard(h)


class FamGenerator(SourceProcess):
    def __init__(self, fam_generator=None, **kwargs):
        super().__init__(**kwargs)
        self.fams_to_process = fam_generator

    def generate_data(self):
        for fam in self.fams_to_process:
            yield fam
        logger.info("all families put to queue")


class HashWorker(BaseProfileBuilderProcess):
    def __init__(self, db_path, **kwargs):
        super().__init__(**kwargs)
        self.db_path = db_path
        self.db = None
        self.hasher = None
        self.handled_queries = None

    def setup(self):
        self.db = Database(self.db_path)
        self.hasher = HogHasher(self.db)
        self.handled_queries = 0

    def handle_input(self, fam):
        logger.info("computing hashes for family {}".format(fam))
        if self.handled_queries > 10000:
            self.db.close()
            logger.info("resetting database handle")
            self.setup()
        logger.info("start chewing on fam {}...".format(fam))
        t0 = time.time()
        hashes = self.hasher.analyze_fam(fam)
        logger.info("... done with fam {}. Took {} sec".format(fam, time.time() - t0))
        self.handled_queries += 1
        return hashes

    def finalize(self):
        self.db.close()


class Collector(BaseProfileBuilderProcess):
    def __init__(self, output_path, **kwargs):
        super().__init__(**kwargs)
        self.output_path = output_path

    def setup(self):
        self.lsh_builder = LSHBuilder(self.output_path, mode="a")

    def handle_input(self, hashes):
        logger.info(
            "build lsh for {} hashes ({})".format(len(hashes), next(iter(hashes)))
        )
        self.lsh_builder.add_minhashes(hashes.items())

    def finalize(self):
        self.lsh_builder.close()


def generator_of_unprocessed_fams(db_path, lsh_path=None):
    def _get_nr_families(db_path):
        db = Database(db_path)
        nr_hogs = db.get_nr_toplevel_hogs()
        db.close()
        logger.info("Found {} families to process".format(nr_hogs))
        return nr_hogs

    def _load_unprocessed_fams(db_path, lsh_path):
        db = Database(db_path)
        with tables.open_file(lsh_path, "r") as lsh_h5:
            processed_fams = set(db.parse_hog_id(x) for x in lsh_h5.get_node("/hogids"))
        remaining = set(range(1, db.get_nr_toplevel_hogs() + 1)) - processed_fams
        db.close()
        logger.info("Found {} unprocessed families".format(len(remaining)))
        return remaining

    if lsh_path is None or not os.path.exists(lsh_path):
        fams_to_process = range(1, _get_nr_families(db_path) + 1)
    else:
        fams_to_process = _load_unprocessed_fams(db_path, lsh_path)
    return fams_to_process


def compute_minhashes_for_db(db_path, output_path, nr_procs=None):
    pipeline = Pipeline()
    if nr_procs is None:
        nr_procs = multiprocessing.cpu_count()

    fams_to_process = generator_of_unprocessed_fams(db_path, output_path)
    pipeline.add_stage(Stage(FamGenerator, nr_procs=1, fam_generator=fams_to_process))
    pipeline.add_stage(Stage(HashWorker, nr_procs=nr_procs, db_path=db_path))
    pipeline.add_stage(Stage(Collector, nr_procs=1, output_path=output_path))
    print("setup pipeline, about to start it.")
    pipeline.run()
    print("finished with computing the MinHashLSH for {}".format(db_path))


def compare_versions(output_file, target_path, *old_path):
    lsh = LSHBuilder(target_path, mode="r")
    lsh.compute_lsh()
    with tables.open_file(
        output_file,
        "w",
        filters=tables.Filters(
            complevel=1, complib="blosc", shuffle=True, fletcher32=True
        ),
    ) as h5_map:
        tab = h5_map.create_table(
            "/",
            "hogmap",
            description=numpy.dtype(
                [("Old", "S255"), ("New", "S255"), ("Jaccard", "f4")]
            ),
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
            for old_id, old_hashvals in tqdm(
                itertools.zip_longest(old.hogids, old.hashes), total=len(old.hogids)
            ):
                minhash = LeanMinHash256(hashvalues=old_hashvals)
                candidates = sorted(lsh.query(old_id, minhash), key=lambda x: -x[2])
                logger.debug("old_id: {}: candidates: {}".format(old_id, candidates))
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


def build_lookup(target_db, *old_dbs, nr_procs=None):
    compute_minhashes_for_db(target_db, target_db + ".hog-lsh.h5", nr_procs=nr_procs)
    for old in old_dbs:
        compute_minhashes_for_db(old, old + "hog-lsh.h5", nr_procs=nr_procs)
