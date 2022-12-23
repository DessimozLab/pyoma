from unittest import TestCase
import numpy
from pyoma.browser import xref_contrib


class GeneEntriesTest(TestCase):
    def test_contiguous(self):
        g = xref_contrib.GeneEntries(numpy.arange(7, 12), 8)
        slices = g.entrynr_slices()
        self.assertEqual(1, len(slices))
        self.assertEqual(7, slices[0].start)
        self.assertEqual(12, slices[0].stop)

    def test_gapped(self):
        g = xref_contrib.GeneEntries(numpy.array([1, 3, 4, 6, 7, 8, 10], dtype="i4"), 4)
        slices = g.entrynr_slices()
        self.assertEqual(4, len(slices))
        self.assertEqual(10, slices[3].start)
        self.assertEqual(11, slices[3].stop)


class SpliceHelperTest(TestCase):
    def test_round_robin(self):
        # entry numbers:   1  2  3  4  5  6  7   8  9, 10, 11  12  13
        alt = numpy.array([0, 2, 2, 0, 6, 6, 6, 10, 9, 10, 11, 10, 11], dtype="i4")
        expect_main = [1, 2, 4, 6, 10, 9, 11]
        expect_sets = [{1}, {2, 3}, {4}, {5, 6, 7}, {8, 10, 12}, {9}, {11, 13}]
        slh = xref_contrib.SpliceVariantHelper(None, alt_array=alt)
        genes = list(slh.iter_genes())
        self.assertEqual(len(expect_sets), len(genes))
        for i in range(len(expect_sets)):
            self.assertEqual(expect_main[i], genes[i].main)
            self.assertEqual(expect_sets[i], set(genes[i].enrs))
