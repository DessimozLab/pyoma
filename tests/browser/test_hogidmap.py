import logging
import random
import unittest

import numpy

from pyoma.browser.db import Database
from pyoma.browser import hogidmap
import tempfile
import os
import math

from .test_db import find_path_to_test_db


class HogIdMapTester(unittest.TestCase):
    db = None

    @classmethod
    def setUpClass(cls):
        path = find_path_to_test_db("TestDb.h5")
        cls.db = Database(path)

    @classmethod
    def tearDownClass(cls):
        cls.db.get_hdf5_handle().close()

    def setUp(self) -> None:
        with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as fh:
            self.hogmapfn = fh.name

    def tearDown(self) -> None:
        try:
            os.remove(self.hogmapfn)
        except OSError:
            pass

    def yield_hog_hashes_for_family_range(self, start=1, stop=40):
        hasher = hogidmap.HogHasher(self.db)
        for fam in range(start, stop):
            yield from hasher.analyze_fam(fam).items()

    def test_query_exact_identical_groups_will_return_same_hogids(self):
        hash_fam_40 = list(self.yield_hog_hashes_for_family_range(stop=40))
        lsh = hogidmap.LSHBuilder(self.hogmapfn, "w")
        lsh.add_minhashes(hash_fam_40)
        lsh.compute_lsh()
        lsh.close()

        lshq = hogidmap.LSHBuilder(self.hogmapfn, "r")
        for hogid, hashes in hash_fam_40:
            best = sorted(lshq.query(hogid, hashes), key=lambda x: -x[2])[0]
            self.assertEqual(
                best[0], best[1], "wrong mapping for {}: {}".format(hogid, best[1])
            )

    def test_stop_and_continue_building_db_works(self):
        hash_fam_40 = list(self.yield_hog_hashes_for_family_range(stop=40))
        lsh = hogidmap.LSHBuilder(self.hogmapfn, "w")
        split = len(hash_fam_40) // 2
        lsh.add_minhashes(hash_fam_40[:split])
        lsh.close()  # close and reopen for append
        lsh = hogidmap.LSHBuilder(self.hogmapfn, "a")
        lsh.add_minhashes(hash_fam_40[split:])
        lsh.compute_lsh()

        for hogid, hashes in hash_fam_40:
            best = sorted(lsh.query(hogid, hashes), key=lambda x: -x[2])[0]
            self.assertEqual(
                best[0], best[1], "wrong mapping for {}: {}".format(hogid, best[1])
            )

    def test_inexisting_hog_has_no_match(self):
        lsh = hogidmap.LSHBuilder(self.hogmapfn, "w")
        lsh.add_minhashes(self.yield_hog_hashes_for_family_range(stop=40))
        for hogid, hashes in self.yield_hog_hashes_for_family_range(start=41, stop=44):
            self.assertEqual(0, len(list(lsh.query(hogid, hashes))))

    def test_largely_overlap_matches(self):
        lsh = hogidmap.LSHBuilder(self.hogmapfn, "w")
        lsh.add_minhashes(self.yield_hog_hashes_for_family_range(start=450, stop=480))
        lsh.compute_lsh()

        query_sub_fam = "HOG:0000474.1d.2f"
        all_member_of_subfam = self.db.member_of_hog_id(query_sub_fam)
        nr_matching = math.ceil(len(all_member_of_subfam) * 0.85)
        logging.warning(
            "building samples of {} ids from pool of {}".format(
                nr_matching, len(all_member_of_subfam)
            )
        )

        no_match = 0
        for it in range(10):
            query_members = numpy.random.choice(
                all_member_of_subfam, nr_matching, replace=False
            )
            minhash = hogidmap.MinHash256()
            for memb in query_members:
                minhash.update(memb["CanonicalId"])
            res = sorted(lsh.query("query", minhash), key=lambda x: -x[2])
            if len(res) == 0:
                no_match += 1
                logging.warning("no match for {}".format(query_members))
                continue
            self.assertEqual(query_sub_fam, res[0][1].decode())
        self.assertLessEqual(no_match, 1)
