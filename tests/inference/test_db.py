import os

import unittest
from pyoma.inference.tabledef import BestMatchesFmt, GenomeFmt
import pyoma.inference.db
import tables
import numpy as np
import tempfile


def get_test_data():
    return (
        np.array(
            [
                (1, 1, 251.32, 22.4, 9.11, 422, True, False, True),
                (1, 5, 2251.32, 2.4, 3.11, 852, True, True, True),
                (1, 6, 151.32, 122.4, 59.11, 223, False, False, False),
                (2, 4, 211.55, 33.1, 10.1, 599, True, True, True),
                (2, 8, 1111.52, 13.1, 10.1, 799, False, False, True),
                (5, 6, 222, 31.7332, 222.1, 554, True, True, True),
            ],
            dtype=tables.dtype_from_descr(BestMatchesFmt),
        ),
        np.array([0, 3, 3, 5, 5, 5, 5, 5, 5, 6]).reshape(5, 2),
    )


class RelationsOfEntryTest(unittest.TestCase):
    def setUp(self):
        data, idx = get_test_data()
        self.indata = data[np.where(data["EntryNr1"] == 1)]

    def test_contains(self):
        all_pairs = pyoma.inference.db.RelationsOfEntry(self.indata)
        self.assertFalse(4 in all_pairs)
        for en2 in (1, 5, 6):
            self.assertTrue(en2 in all_pairs, str(en2) + " is not in all_pairs")

        spairs = pyoma.inference.db.StablePairsOfEntry(self.indata)
        self.assertFalse(6 in spairs)
        for en2 in (1, 5):
            self.assertTrue(en2 in spairs)

        vpairs = pyoma.inference.db.VPairsOfEntry(self.indata)
        for en2 in (2, 1, 6):
            self.assertFalse(en2 in vpairs)
        self.assertTrue(5 in vpairs)

    def test_iterator(self):
        spairs = pyoma.inference.db.StablePairsOfEntry(self.indata)
        self.assertEqual(
            [x["EntryNr2"] for x in spairs], [1, 5], "Unexpected members in SPairs"
        )

    def test_matches(self):
        vpairs = pyoma.inference.db.VPairsOfEntry(self.indata)
        self.assertEqual(
            1, len(vpairs.relations()), "Unexpected number of VPs returned"
        )


class GenomePairTest(unittest.TestCase):
    def setUp(self):
        matches, itab = get_test_data()
        itab = np.array([0, 3, 3, 5, 5, 5, 5, 5, 5, 6]).reshape(5, 2)
        self.genome_pair = pyoma.inference.db.GenomePair(matches, itab)

    def test_get_entry(self):
        prot1_matches = [
            (m["EntryNr1"], m["EntryNr2"]) for m in self.genome_pair[1]["ALL"]
        ]
        self.assertEqual(
            [(1, 1), (1, 5), (1, 6)],
            prot1_matches,
            "unexpected matches with ALL selector",
        )
        prot2_matches = [
            (m["EntryNr1"], m["EntryNr2"]) for m in self.genome_pair[2]["SP"]
        ]
        self.assertEqual([(2, 4)], prot2_matches, "unexpected matches with SP selector")
        self.assertEqual(
            0,
            len(self.genome_pair[3]["ALL"].relations()),
            "Unexpected matches, should be empty",
        )

    def test_keep_different_reltypes_in_sync(self):
        prot1_before = [z["EntryNr2"] for z in self.genome_pair[1]["BM"]]
        self.assertEqual([1, 5], prot1_before, "test-setup failed")
        self.genome_pair[1]["VP"] = [1, 6]
        prot2_after = [z["EntryNr2"] for z in self.genome_pair[1]["BM"]]
        self.assertEqual([1, 5, 6], prot2_after, "change of VP should set BM flag")

    def test_fail_if_not_member(self):
        self.assertRaises(
            ValueError, self.genome_pair[2]["BM"].set_relations, [4, 6, 8]
        )

    def test_iterators_in_sync(self):
        # create iterator before any changes on data
        bm_it = iter(self.genome_pair[1]["BM"])
        # set new bm pairs
        new_matches = [5, 6]
        self.genome_pair[1]["BM"] = new_matches
        matches_from_iterator = [z["EntryNr2"] for z in bm_it]
        self.assertEqual(matches_from_iterator, new_matches)


class OmaDBTest(unittest.TestCase):
    def setUp(self):
        matches, itab = get_test_data()
        genome_summary = np.array(
            [
                (1234, "SE001", 5, 21122, 0, "Synthetic Genome 1", "Test"),
                (1235, "SE002", 8, 15222, 5, "Synthetic Genome 2", "Test"),
            ],
            tables.dtype_from_descr(GenomeFmt),
        )
        self.db_filename = tempfile.mktemp(".h5")
        fd = tables.open_file(self.db_filename, "w")
        se1vs2 = fd.create_group("/Matches/SE001", "SE002", createparents=True)
        fd.create_table(se1vs2, "Relations", BestMatchesFmt, obj=matches)
        fd.create_carray(se1vs2, "ProteinIndex", obj=itab)
        fd.create_table("/", "GenomeSummary", GenomeFmt, obj=genome_summary)
        fd.close()

    def tearDown(self):
        try:
            os.remove(self.db_filename)
        except IOError:
            pass

    def test_without_write_sps(self):
        db = pyoma.inference.db.OmaDB(self.db_filename, mode="a")
        pair = db.matches("SE001", "SE002")
        sp = [z["EntryNr2"] for z in pair[1]["SP"]]
        sp.append(6)
        pair[1]["SP"] = sp
        sp_after = [z["EntryNr2"] for z in pair[1]["SP"]]
        self.assertEqual(sp, sp_after)
        db.close()

    def test_write_vps(self):
        db = pyoma.inference.db.OmaDB(self.db_filename, mode="a")
        vp = [z["EntryNr2"] for z in db.matches("SE001", "SE002")[1]["VP"]]
        vp.append(6)
        db.matches("SE001", "SE002")[1]["VP"] = vp
        db.matches("SE001", "SE002").flush()
        db.close()

        db = pyoma.inference.db.OmaDB(self.db_filename, mode="r")
        vp = [z["EntryNr2"] for z in db.matches("SE001", "SE002")[1]["VP"]]
        self.assertIn(6, vp)
        db.close()
