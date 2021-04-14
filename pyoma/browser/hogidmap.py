import collections
import multiprocessing
import re
import os
from functools import partial

import tables
from datasketch import MinHash, MinHashLSH

from .db import Database
from .hogprofile.build import BaseProfileBuilderProcess, SourceProcess, Pipeline, Stage

MinHash256 = partial(MinHash, num_perm=256)


class HogHasher(object):
    def __init__(self, db: Database):
        self.db = db

    def analyze_fam(self, fam_nr):
        members = self.db.member_of_fam(fam_nr)
        minhashes = collections.defaultdict(MinHash256)
        for e in members:
            hog_id = e["OmaHOG"]
            for p in re.finditer(br"\.", hog_id):
                minhashes[hog_id[: p.start()]].update(e["CanonicalId"])
            minhashes[hog_id].update(e["CanonicalId"])
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
        self._row = len(self.hashes)
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
        self.h5.close()

    def add_minhashes(self, it):
        for key, minhash in it:
            self.hashes.append([minhash.digest()])
            self.lsh.insert(key, minhash)
            self.hogid2row[key] = self._row
            self._row += 1

    def query(self, key, minhash):
        """query with a minhash, returns a list of tuples of
        (query-key, target-key, jaccard)"""
        candidates = self.lsh.query(minhash)
        for c in candidates:
            hashvals = self.hashes[self.hogid2row[c]]
            h = MinHash(seed=1, hashvalues=hashvals)
            yield key, c, minhash.jaccard(h)


class FamGenerator(SourceProcess):
    def __init__(self, db_path, **kwargs):
        super().__init__(**kwargs)
        self.db_path = db_path

    def generate_data(self):
        db = Database(self.db_path)
        nr_hogs = db.get_nr_toplevel_hogs()
        db.close()
        for fam in range(1, nr_hogs + 1):
            yield fam
        print("all families put to queue")


class HashWorker(BaseProfileBuilderProcess):
    def __init__(self, db_path, **kwargs):
        super().__init__(**kwargs)
        self.db_path = db_path

    def setup(self):
        self.db = Database(self.db_path)
        self.hasher = HogHasher(self.db)

    def handle_input(self, fam):
        return self.hasher.analyze_fam(fam)

    def finalize(self):
        self.db.close()


class Collector(BaseProfileBuilderProcess):
    def __init__(self, output_path, **kwargs):
        super().__init__(**kwargs)
        self.output_path = output_path

    def setup(self):
        self.lsh_builder = LSHBuilder(self.output_path, mode="w")

    def handle_input(self, hashes):
        self.lsh_builder.add_minhashes(hashes.items())

    def finalize(self):
        self.lsh_builder.close()


def compute_minhashes_for_db(db_path, output_path, nr_procs=None):
    pipeline = Pipeline()
    if nr_procs is None:
        nr_procs = multiprocessing.cpu_count()

    pipeline.add_stage(Stage(FamGenerator, nr_procs=1, db_path=db_path))
    pipeline.add_stage(Stage(HashWorker, nr_procs=nr_procs, db_path=db_path))
    pipeline.add_stage(Stage(Collector, nr_procs=1, output_path=output_path))
    print("setup pipeline, about to start it.")
    pipeline.run()
    print("finished with computing the MinHashLSH for {}".format(db_path))


def build_lookup(target_db, *old_dbs):
    compute_minhashes_for_db(target_db, target_db + ".hog-lsh.h5")
    for old in old_dbs:
        compute_minhashes_for_db(old, old + "hog-lsh.h5")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="find related hogs in new version")
    parser.add_argument("--target", required=True, help="Path to the target database")
    parser.add_argument(
        "--old", nargs="+", help="Path to databases that should be mapped"
    )
    conf = parser.parse_args()
    print(conf)

    build_lookup(conf.target, conf.old)
    print("done... bye bye")
