import collections
import re
import os
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
