import os
import tempfile
import unittest

import numpy
import tables

from pyoma.browser.idmapper import GeneNamesLookup


def _build_h5(path, entries):
    """Write a minimal XRefIndex/GeneNames + GeneNames_lookup HDF5 at *path*.

    *entries* is a list of ``(name_bytes, [(entry_nr, xref_row), ...])`` pairs.
    Names must already be lowercase.
    """
    entries = sorted(entries, key=lambda x: x[0])  # sort by name, as the real builder does
    max_key_len = max(len(name) for name, _ in entries)
    key_dtype = numpy.dtype([("XRefId", f"S{max_key_len}"), ("Offset", "i4"), ("Length", "i4")])
    lookup_dtype = numpy.dtype([("EntryNr", "i4"), ("XRefRow", "i4")])

    with tables.open_file(path, "w") as h5:
        h5.create_group("/", "XRefIndex")
        key_tab = h5.create_table("/XRefIndex", "GeneNames", description=key_dtype)
        lookup_tab = h5.create_table("/XRefIndex", "GeneNames_lookup", description=lookup_dtype)

        off = 0
        for name, proteins in entries:
            for enr, xrr in proteins:
                lookup_tab.append([(enr, xrr)])
            key_tab.append([(name, off, len(proteins))])
            off += len(proteins)

        key_tab.flush()
        lookup_tab.flush()


class GeneNamesLookupTest(unittest.TestCase):
    # Test data: (lowercase_name, [(entry_nr, xref_row), ...])
    ENTRIES = [
        (b"brca1", [(1, 10), (5, 50)]),
        (b"brca2", [(2, 20)]),
        (b"myc", [(3, 30)]),
        (b"mycn", [(4, 40)]),
        (b"tp53", [(6, 60)]),
    ]

    def setUp(self):
        fd, self.h5_path = tempfile.mkstemp(suffix=".h5")
        os.close(fd)
        _build_h5(self.h5_path, self.ENTRIES)
        self.h5 = tables.open_file(self.h5_path, "r")
        self.lookup = GeneNamesLookup(self.h5)

    def tearDown(self):
        self.h5.close()
        os.unlink(self.h5_path)

    # --- __contains__ ---

    def test_contains_exact(self):
        self.assertIn("brca1", self.lookup)

    def test_contains_case_insensitive(self):
        self.assertIn("BRCA1", self.lookup)
        self.assertIn("Tp53", self.lookup)

    def test_not_contains(self):
        self.assertNotIn("brca3", self.lookup)

    # --- count ---

    def test_count_multi(self):
        self.assertEqual(2, self.lookup.count("brca1"))

    def test_count_single(self):
        self.assertEqual(1, self.lookup.count("tp53"))

    def test_count_absent(self):
        self.assertEqual(0, self.lookup.count("notexist"))

    def test_count_case_insensitive(self):
        self.assertEqual(2, self.lookup.count("BRCA1"))

    # --- get_matching_xref_row_nrs ---

    def test_get_xref_rows(self):
        rows = self.lookup.get_matching_xref_row_nrs("brca1")
        self.assertEqual(set(rows), {10, 50})

    def test_get_xref_rows_with_entrynr_range(self):
        rows = self.lookup.get_matching_xref_row_nrs("brca1", entrynr_range=(1, 5))
        self.assertEqual(list(rows), [10])

    def test_get_xref_rows_missing_raises(self):
        with self.assertRaises(KeyError):
            self.lookup.get_matching_xref_row_nrs("notexist")

    # --- starts_with ---

    def test_starts_with_prefix(self):
        results = self.lookup.starts_with("brca")
        self.assertEqual([name for name, _, _ in results], ["brca1", "brca2"])

    def test_starts_with_case_insensitive(self):
        results = self.lookup.starts_with("BRCA")
        self.assertEqual([name for name, _, _ in results], ["brca1", "brca2"])

    def test_starts_with_includes_count(self):
        results = {name: count for name, count, _ in self.lookup.starts_with("brca")}
        self.assertEqual(results["brca1"], 2)
        self.assertEqual(results["brca2"], 1)

    def test_starts_with_xref_row(self):
        results = {name: xrr for name, _, xrr in self.lookup.starts_with("myc")}
        # xref_row is the first lookup entry for each name
        self.assertEqual(results["myc"], 30)
        self.assertEqual(results["mycn"], 40)

    def test_starts_with_no_match(self):
        self.assertEqual([], self.lookup.starts_with("xyz"))

    def test_starts_with_full_name(self):
        results = self.lookup.starts_with("tp53")
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0][0], "tp53")

    # --- contains ---

    def test_contains_substring(self):
        names = [name for name, _, _ in self.lookup.contains("rc")]
        self.assertEqual(names, ["brca1", "brca2"])

    def test_contains_case_insensitive(self):
        names = [name for name, _, _ in self.lookup.contains("RC")]
        self.assertEqual(names, ["brca1", "brca2"])

    def test_contains_no_match(self):
        self.assertEqual([], self.lookup.contains("xyz"))

    def test_contains_includes_count_and_xref_row(self):
        results = {name: (count, xrr) for name, count, xrr in self.lookup.contains("53")}
        self.assertIn("tp53", results)
        self.assertEqual(results["tp53"], (1, 60))


if __name__ == "__main__":
    unittest.main()
