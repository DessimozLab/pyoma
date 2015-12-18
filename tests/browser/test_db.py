import unittest
import numpy
from pyoma.browser.db import *
from pyoma.browser import tablefmt


class TestHelperFunctions(unittest.TestCase):
    def test_counter(self):
        self.assertEqual(0, count_elements([]))
        self.assertEqual(3, count_elements('abc'))
        recarray = numpy.zeros(2, dtype=[('A','i4'),('B','f8')])
        self.assertEqual(2, count_elements(recarray))


class DatabaseTests(unittest.TestCase):
    db = None

    @classmethod
    def setUpClass(cls):
        path = "/pub/projects/cbrg-oma-browser/Test.Jul2014/data/OmaServer.h5"
        cls.db = Database(path)

    @classmethod
    def tearDownClass(cls):
        cls.db.get_hdf5_handle().close()

    def test_get_vpairs_of_entry_with_orthologs(self):
        for entry_nr, exp_vps_cnt in [(12, 3), (1, 0), (4,1)]:
            vps = self.db.get_vpairs(entry_nr)
            self.assertTrue(isinstance(vps, numpy.ndarray))
            self.assertEqual(exp_vps_cnt, len(vps))
            self.assertEqual(['EntryNr1', 'EntryNr2', 'RelType'],
                             sorted(vps.dtype.fields.keys()))

    def test_neighborhood_close_to_boundary(self):
        query, window = 3, 7
        neighbors, idx = self.db.neighbour_genes(query, window)
        self.assertEqual(query - 1, idx)
        self.assertEqual(neighbors['EntryNr'][idx], query)
        expected_entry_nrs = numpy.arange(1, query + window + 1, dtype='i4')
        self.assertTrue(numpy.array_equal(expected_entry_nrs,
                                          neighbors['EntryNr']))

    def test_hog_family(self):
        entry = numpy.zeros(1, dtype=tables.dtype_from_descr(tablefmt.ProteinTable))
        entry['OmaHOG'] = b""
        with self.assertRaises(Singleton):
            self.db.hog_family(entry)
        entry['OmaHOG'] = b"HOG:000523"
        self.assertEqual(523, self.db.hog_family(entry))

    def test_orthoxml(self):
        xml = self.db.get_orthoxml(33).decode()
        # we simply check that orthoxml starts with <?xml and ends with an orthoxml tag
        self.assertTrue(xml.startswith('<?xml '))
        self.assertTrue(xml.endswith('</orthoXML>\n'))

    def test_hog_lex_range(self):
        cases = [(b'HOG:001', (b'HOG:001', b'HOG:002')),
                 (b'HOG:001.1a.2b', (b'HOG:001.1a.2b', b'HOG:001.1a.2c'))]
        for hog, rng in cases:
            self.assertEqual(rng, self.db._hog_lex_range(hog))

    def test_hog_members(self):
        cases = [('Eukaryota', 5), ('Fungi', 3), ('Taphrinomycotina',0)]
        for level, exp_member_cnt in cases:
            if exp_member_cnt == 0:
                with self.assertRaises(ValueError):
                    self.db.hog_members(12, level)
            else:
                self.assertEqual(len(self.db.hog_members(12, level)), exp_member_cnt)