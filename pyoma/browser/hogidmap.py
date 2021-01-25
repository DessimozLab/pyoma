import collections
import re
from functools import partial

import tables
from datasketch import MinHash, MinHashLSH

from .db import Database

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
        return minhashes.items()


class LSHBuilder(object):
    def __init__(self, threshold=0.8):
        self.lsh = MinHashLSH(threshold=threshold, num_perm=256)
        self.h5 = self.init_hash_table_file()
        self.hashes = self.h5.get_node("/hashes")
        self.hogid2row = {}
        self._row = 0

    def init_hash_table_file(self):
        h5 = tables.open_file(
            "hogmap-hashes.h5",
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
        for key, minhash in it:
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
