import collections
import re
from functools import partial
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import threading
import os
import tables
from datasketch import MinHash, MinHashLSH
from datasketch.experimental.aio.lsh import AsyncMinHashLSH
from .db import Database

MinHash256 = partial(MinHash, num_perm=256)


class HogHasher(object):
    def __init__(self, dbfn: str):
        print("init hoghasher: {}".format(self))
        self.dbfn = dbfn

    def __call__(self, fam_nr, **kwargs):
        try:
            db = self.db
        except AttributeError:
            self.db = Database(self.dbfn)
            db = self.db

        print("fam {}".format(fam_nr))
        members = db.member_of_fam(fam_nr)
        minhashes = collections.defaultdict(MinHash256)
        for e in members:
            hog_id = e["OmaHOG"]
            for p in re.finditer(br"\.", hog_id):
                minhashes[hog_id[: p.start()]].update(e["CanonicalId"])
            minhashes[hog_id].update(e["CanonicalId"])
        return minhashes.items()


def doubler(value):
    return ("hog", MinHash256(set(value))), ("hog2", MinHash256(set(value * 3)))


class LSHBuilder(object):
    def __init__(self, threshold=0.8, hash_fname=None):
        self.lsh = MinHashLSH(threshold=threshold, num_perm=256)
        self.h5 = self.init_hash_table_file(hash_fname)
        self.hashes = self.h5.get_node("/hashes")
        self.hogid2row = {}
        self._row = 0

    def init_hash_table_file(self, hash_fname=None):
        if hash_fname is None:
            hash_fname = os.path.join(
                os.getenv("DARWIN_BROWSER_SCRATCH_PATH", ""), "hogmap-hashes.h5"
            )
        h5 = tables.open_file(
            hash_fname,
            "w",
            filters=tables.Filters(
                complevel=1, complib="blosc", shuffle=True, fletcher32=True
            ),
        )
        h5.create_earray(
            "/", "hashes", atom=tables.Int64Atom(), shape=(0, 256), expectedrows=1e6
        )
        h5.create_vlarray("/", "lsh_obj", atom=tables.ObjectAtom())
        return h5

    def close(self):
        lsh_obj_arr = self.h5.get_node("/lsh_obj")
        lsh_obj_arr.append(self.lsh)
        lsh_obj_arr.append(self.hogid2row)
        self.hashes.flush()
        self.h5.close()

    def add_minhashes(self, it):
        mh = []
        for row, (key, minhash) in enumerate(it, start=self._row):
            self.lsh.insert(key, minhash)
            self.hogid2row[key] = row
            mh.append(minhash)
        self.hashes.append(mh)

    def add_minhash(self, key, minhash):
        self.lsh.insert(key, minhash)
        self.hashes.append([minhash.digest()])
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


def build_lsh(dbfn):
    print("init lsh...")
    lsh = LSHBuilder()
    print("done init lsh")
    db = Database(dbfn)
    nr_hogs = db.get_nr_toplevel_hogs()
    db.close()
    _run(lsh, dbfn, nr_hogs)


def _run(lsh, dbfn, nr_hogs):
    worker = HogHasher(dbfn)
    families = range(1, nr_hogs + 1)
    with ProcessPoolExecutor(max_workers=4) as exec:
        for result in exec.map(doubler, families):
            print(result)
            # lsh.add_minhash(*result)
    return lsh
