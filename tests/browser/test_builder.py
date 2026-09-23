from __future__ import division, print_function, unicode_literals

import gzip
import json
import os
import random
import shutil
import tempfile
import unittest

import numpy
import numpy.testing
import pandas
import tables

from pyoma.browser import tablefmt
from pyoma.browser.KmerEncoder import KmerEncoder, DIGITS_AA
from pyoma.browser.exceptions import DBConsistencyError
from pyoma.browser.build.builder import (
    DBBuilder,
    OmaGroupsProvider,
    XrefStorer,
    identify_close_paralogs,
    load_homoeologs_from_tsv,
    log_and_replace_error_handler,
)


def make_taxonomy_dtype_row(**kwargs):
    row = numpy.zeros(1, dtype=tables.dtype_from_descr(tablefmt.TaxonomyTable))[0]
    for k, v in kwargs.items():
        row[k] = v
    return row


class LogAndReplaceErrorHandlerTest(unittest.TestCase):
    def test_replaces_bad_byte_with_replacement_char(self):
        try:
            b"\xff".decode("utf-8")
        except UnicodeDecodeError as e:
            repl, end = log_and_replace_error_handler(e)
            self.assertEqual("�", repl)
            self.assertEqual(e.end, end)
        else:
            self.fail("expected UnicodeDecodeError")


class OmaGroupsProviderTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.source = os.path.join(self.tmpdir, "groups.json")
        with open(self.source, "w") as fh:
            json.dump({"HUMAN": {"1": 42, "2": 7}}, fh)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_known_genome_and_nr(self):
        prov = OmaGroupsProvider(self.source)
        self.assertEqual(42, prov.get_oma_group("HUMAN", 1))
        self.assertEqual(7, prov.get_oma_group("HUMAN", 2))

    def test_unknown_nr_returns_zero(self):
        prov = OmaGroupsProvider(self.source)
        self.assertEqual(0, prov.get_oma_group("HUMAN", 999))

    def test_unknown_genome_returns_zero(self):
        prov = OmaGroupsProvider(self.source)
        self.assertEqual(0, prov.get_oma_group("MOUSE", 1))

    def test_nr_as_int_falls_back_correctly(self):
        # json only has string keys; the first lookup with an int key fails
        # with KeyError, and get_oma_group falls back to trying the same
        # int key again (which also fails), so it must return 0 -- but if
        # some data was written with proper int-vs-str duplication this
        # branch would be exercised. Here we mainly confirm it does not
        # raise and returns 0 for a lookup that only exists as str key
        # when queried with a mismatched type via nr already being int
        # (the normal call path from add_proteins passes an int).
        prov = OmaGroupsProvider(self.source)
        self.assertEqual(42, prov.get_oma_group("HUMAN", 1))

    def test_none_source_raises_typeerror_on_lookup(self):
        # BUG-PIN: OmaGroupsProvider(None) sets self.data = None. get_oma_group
        # then does `self.data[str(genome)]` which raises TypeError (not
        # subscriptable), which is *not* caught by the `except KeyError`
        # clause in builder.py OmaGroupsProvider.get_oma_group. This is
        # called as OmaGroupsProvider(conf.oma_groups) in build/main.py,
        # where conf.oma_groups may be None/unset -- pinning current
        # (likely buggy) behavior rather than fixing it here.
        prov = OmaGroupsProvider(None)
        self.assertIsNone(prov.data)
        with self.assertRaises(TypeError):
            prov.get_oma_group("HUMAN", 1)


class IdentifyClseParalogsTest(unittest.TestCase):
    def _make_df(self):
        # EntryNr2 == 100 is shared by three EntryNr1 (close paralogs)
        # EntryNr2 == 101 is shared by two EntryNr1
        # EntryNr2 == 102 is a singleton (no close paralogs)
        rows = [
            (1, 100),
            (2, 100),
            (3, 100),
            (4, 101),
            (5, 101),
            (6, 102),
        ]
        return pandas.DataFrame(rows, columns=["EntryNr1", "EntryNr2"])

    def _expected_pairs(self):
        return {(1, 2), (1, 3), (2, 3), (4, 5)}

    def test_join_based_method(self):
        df = self._make_df()
        res = identify_close_paralogs(df, join_threshold_mb=500)
        pairs = set(zip(res["EntryNr1"], res["EntryNr2"]))
        self.assertEqual(self._expected_pairs(), pairs)
        self.assertTrue((res["EntryNr1"] < res["EntryNr2"]).all())
        self.assertEqual(len(pairs), len(res))
        self.assertTrue(list(zip(res["EntryNr1"], res["EntryNr2"])) == sorted(zip(res["EntryNr1"], res["EntryNr2"])))

    def test_groupby_based_method(self):
        df = self._make_df()
        res = identify_close_paralogs(df, join_threshold_mb=0)
        pairs = set(zip(res["EntryNr1"], res["EntryNr2"]))
        self.assertEqual(self._expected_pairs(), pairs)
        self.assertTrue((res["EntryNr1"] < res["EntryNr2"]).all())
        self.assertEqual(len(pairs), len(res))
        self.assertTrue(list(zip(res["EntryNr1"], res["EntryNr2"])) == sorted(zip(res["EntryNr1"], res["EntryNr2"])))

    def test_both_methods_agree(self):
        df = self._make_df()
        join_res = identify_close_paralogs(df, join_threshold_mb=500)
        groupby_res = identify_close_paralogs(df, join_threshold_mb=0)
        self.assertEqual(
            set(zip(join_res["EntryNr1"], join_res["EntryNr2"])),
            set(zip(groupby_res["EntryNr1"], groupby_res["EntryNr2"])),
        )

    def test_no_shared_entrynr2_groupby_method_gives_empty_result(self):
        # the groupby-based method handles the "no close paralogs at all"
        # case correctly. Note: when there are truly no shared EntryNr2
        # values, the estimated memory of the join-based method is always 0
        # (est_pairs=0), so `join_threshold_mb=0` is not enough to force the
        # groupby branch (0 <= 0 still picks the join method) -- use a
        # negative threshold to force it.
        df = pandas.DataFrame({"EntryNr1": [1, 2, 3], "EntryNr2": [10, 20, 30]})
        res_groupby = identify_close_paralogs(df, join_threshold_mb=-1)
        self.assertEqual(0, len(res_groupby))

    def test_no_shared_entrynr2_join_method_raises_BUG(self):
        # BUG-PIN: pyoma/browser/build/builder.py::identify_close_paralogs
        # (join-based branch, ~lines 374-377). When there are no shared
        # EntryNr2 values at all, the intermediate `cp` DataFrame after the
        # `cp[cp["EntryNr1"] < cp["EntryNr1_2"]]` filter is empty, and
        # `.drop_duplicates(ignore_index=True)` on an *empty* DataFrame does
        # not actually clear the (inherited) index name "EntryNr2" in this
        # pandas version -- even though `ignore_index=True` is supposed to
        # give a fresh default index. The subsequent
        # `.rename(columns={"EntryNr1_2": "EntryNr2"})` then produces a frame
        # with BOTH a column and an index level named "EntryNr2", which makes
        # the final `cp.sort_values(by=["EntryNr1", "EntryNr2"], ...)` raise
        # `ValueError: 'EntryNr2' is both an index level and a column label,
        # which is ambiguous.` This only manifests for the empty-result edge
        # case (e.g. a genome pair with zero close paralogs) via the
        # join-based method (`join_threshold_mb` large enough, the default
        # in production use). Pinning current behavior; not fixed here.
        df = pandas.DataFrame({"EntryNr1": [1, 2, 3], "EntryNr2": [10, 20, 30]})
        with self.assertRaises(ValueError):
            identify_close_paralogs(df, join_threshold_mb=500)


class BufferedTableWriterTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.dtype = numpy.dtype([("A", "i4"), ("B", "f8")])

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _new_db_and_table(self, name="t.h5", tabname="T"):
        path = os.path.join(self.tmpdir, name)
        db = DBBuilder(path)
        db.__enter__()
        self.addCleanup(db.__exit__, None, None, None)
        tab = db.h5.create_table("/", tabname, description=self.dtype)
        return db, tab

    def test_add_one_and_add_list_of_tuples(self):
        from pyoma.browser.build.builder import BufferedTableWriter

        db, tab = self._new_db_and_table()
        w = BufferedTableWriter(tab, dtype=self.dtype)
        w.add_one((1, 1.5))
        w.add([(2, 2.5), (3, 3.5)])
        w.flush()
        data = tab.read()
        self.assertEqual(3, len(data))
        numpy.testing.assert_array_equal([1, 2, 3], data["A"])
        numpy.testing.assert_allclose([1.5, 2.5, 3.5], data["B"])
        self.assertEqual(3, w.total_written)

    def test_add_pandas_dataframe(self):
        from pyoma.browser.build.builder import BufferedTableWriter

        db, tab = self._new_db_and_table()
        w = BufferedTableWriter(tab, dtype=self.dtype)
        df = pandas.DataFrame({"A": [5, 6], "B": [1.1, 2.2]})
        w.add(df)
        w.flush()
        data = tab.read()
        self.assertEqual(2, len(data))
        numpy.testing.assert_array_equal([5, 6], sorted(data["A"]))

    def test_add_numpy_structured_array(self):
        from pyoma.browser.build.builder import BufferedTableWriter

        db, tab = self._new_db_and_table()
        w = BufferedTableWriter(tab, dtype=self.dtype)
        arr = numpy.array([(7, 7.7), (8, 8.8)], dtype=self.dtype)
        w.add(arr)
        w.flush()
        data = tab.read()
        self.assertEqual(2, len(data))
        numpy.testing.assert_array_equal([7, 8], data["A"])

    def test_sort_key_sorts_on_flush(self):
        from pyoma.browser.build.builder import BufferedTableWriter

        db, tab = self._new_db_and_table()
        w = BufferedTableWriter(tab, sort_key=["A"], dtype=self.dtype)
        w.add([(3, 0.0), (1, 0.0), (2, 0.0)])
        w.flush()
        data = tab.read()
        numpy.testing.assert_array_equal([1, 2, 3], data["A"])

    def test_autoflush_on_buffer_size(self):
        from pyoma.browser.build.builder import BufferedTableWriter

        db, tab = self._new_db_and_table()
        w = BufferedTableWriter(tab, buffer_size=3, dtype=self.dtype)
        w.add([(1, 0.0), (2, 0.0), (3, 0.0)])
        # buffer_cnt reached buffer_size -> auto flush already happened
        self.assertEqual(0, w.buffer_cnt)
        self.assertEqual(3, w.total_written)
        self.assertEqual(3, len(tab))

    def test_add_empty_list_is_noop(self):
        from pyoma.browser.build.builder import BufferedTableWriter

        db, tab = self._new_db_and_table()
        w = BufferedTableWriter(tab, dtype=self.dtype)
        w.add([])
        self.assertEqual(0, w.buffer_cnt)
        w.flush()
        self.assertEqual(0, len(tab))


class XrefStorerTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_basic_single_file_flow(self):
        path = os.path.join(self.tmpdir, "xref.h5")
        with XrefStorer(path, mode="w") as st:
            st.add_source_xref(1, "ID1", "id")
            st.add_source_xref(1, "AC1", "ac")
            st.add_ec(1, "EC1")

        with tables.open_file(path, mode="r") as h5:
            xref = h5.root.XRef
            self.assertEqual(2, len(xref))
            source_enum = xref.get_enum("XRefSource")
            rows = {(r["EntryNr"], r["XRefId"].decode(), r["XRefSource"]) for r in xref.iterrows()}
            self.assertIn((1, "ID1", source_enum["SourceID"]), rows)
            self.assertIn((1, "AC1", source_enum["SourceAC"]), rows)

            ec = h5.root.Annotations.EC
            self.assertEqual(1, len(ec))
            self.assertEqual(1, ec[0]["EntryNr"])
            self.assertEqual(b"EC1", ec[0]["ECacc"])

    def test_rows_come_out_sorted(self):
        path = os.path.join(self.tmpdir, "xref_sorted.h5")
        with XrefStorer(path, mode="w") as st:
            # add out of (EntryNr, XRefSource, XRefId, Verification) order
            st.add_xref(3, 10, "b", 0, 1.0)
            st.add_xref(1, 20, "a", 0, 1.0)
            st.add_xref(1, 10, "z", 0, 1.0)
            st.add_xref(2, 5, "m", 0, 1.0)

        with tables.open_file(path, mode="r") as h5:
            xref = h5.root.XRef
            keys = [(r["EntryNr"], r["XRefSource"], r["XRefId"], r["Verification"]) for r in xref.iterrows()]
            self.assertEqual(sorted(keys), keys)

    def test_multi_files_rotation(self):
        base = os.path.join(self.tmpdir, "multi.h5")
        # NOTE: rotation is only checked against *flushed* row counts
        # (XrefFileWriter.total_rows sums BufferedTableWriter.total_written,
        # which only advances on flush()). With the default buffer_size
        # (500_000) a handful of added rows would never trigger a flush and
        # therefore never rotate, even past max_rows_per_file. Use
        # buffer_size=1 here so every add() flushes immediately and rotation
        # can actually be observed at this tiny scale.
        with XrefStorer(base, mode="w", multi_files=True, max_rows_per_file=2, buffer_size=1) as st:
            for i in range(1, 6):
                st.add_source_xref(i, f"ID{i}", "id")

        # 5 rows at max_rows_per_file=2 rotate into 3 files (2, 2, 1 rows)
        f0 = os.path.join(self.tmpdir, "multi_000.h5")
        f1 = os.path.join(self.tmpdir, "multi_001.h5")
        f2 = os.path.join(self.tmpdir, "multi_002.h5")
        self.assertTrue(os.path.exists(f0))
        self.assertTrue(os.path.exists(f1))
        self.assertTrue(os.path.exists(f2))
        self.assertFalse(os.path.exists(os.path.join(self.tmpdir, "multi_003.h5")))

        total = 0
        for fn in (f0, f1, f2):
            with tables.open_file(fn, mode="r") as h5:
                total += len(h5.root.XRef)
        self.assertEqual(5, total)

    def test_invalid_index_col_raises(self):
        path = os.path.join(self.tmpdir, "bad.h5")
        with self.assertRaises(ValueError):
            XrefStorer(path, index_cols=["NotAColumn"])


class LoadHomoeologsFromTsvTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _genome_row(self, code="HUMA", off=100):
        row = numpy.zeros(1, dtype=tables.dtype_from_descr(tablefmt.GenomeTable))[0]
        row["UniProtSpeciesCode"] = code.encode("utf-8")
        row["EntryOff"] = off
        return row

    def test_basedir_none_returns_empty(self):
        res = load_homoeologs_from_tsv(self._genome_row(), basedir=None)
        self.assertEqual(0, len(res))
        self.assertEqual(tables.dtype_from_descr(tablefmt.PairwiseRelationTable), res.dtype)

    def test_missing_file_returns_empty(self):
        res = load_homoeologs_from_tsv(self._genome_row(code="XXXX"), basedir=self.tmpdir)
        self.assertEqual(0, len(res))

    def test_valid_file_loads_with_offset(self):
        code = "HUMA"
        genome = self._genome_row(code=code, off=100)
        fn = os.path.join(self.tmpdir, f"{code}.tsv.gz")
        with gzip.open(fn, mode="wt") as fh:
            fh.write("1\t1\t9900\t1:1\t0.9\t2.5\n")
            fh.write("2\t2\t9800\tn:1\t0.8\t3.1\n")
        res = load_homoeologs_from_tsv(genome, basedir=self.tmpdir)
        self.assertEqual(2, len(res))
        numpy.testing.assert_array_equal([101, 102], res["EntryNr1"])
        numpy.testing.assert_array_equal([101, 102], res["EntryNr2"])


class FamIndexHelper(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _build_db_with_entries(self, hog_ids):
        path = os.path.join(self.tmpdir, "famidx_{}.h5".format(len(hog_ids)))
        db = DBBuilder(path)
        db.__enter__()
        self.addCleanup(db.__exit__, None, None, None)
        n = len(hog_ids)
        entries = numpy.zeros(n, dtype=tables.dtype_from_descr(tablefmt.ProteinTable))
        entries["EntryNr"] = numpy.arange(1, n + 1)
        db.h5.create_group("/Protein", "Root", createparents=True)
        db.h5.create_table("/Protein", "Entries", obj=entries)
        return db

    def test_no_rel_char_format(self):
        hog_ids = numpy.array(
            [
                b"HOG:0000001.1a",
                b"HOG:0000001.1b",
                b"HOG:0000002",
                b"",
            ],
            dtype="S255",
        )
        db = self._build_db_with_entries(hog_ids)
        db.add_protein_hog_ids(hog_ids)

        entries = db.h5.get_node("/Protein/Entries").read()
        numpy.testing.assert_array_equal(hog_ids, entries["OmaHOG"])

        lookup = db.h5.get_node("/Protein/FamIndex/Lookup").read()
        total = sum(int(row["count"]) for row in lookup)
        self.assertEqual(len(hog_ids), total)

        entry_idx = db.h5.get_node("/Protein/FamIndex/EntryIdx").read()
        self.assertEqual(len(hog_ids), len(entry_idx))
        reordered = hog_ids[entry_idx]

        # family 0 (empty OmaHOG) bucket must be present with the right count
        fam0_count = sum(1 for h in hog_ids if h == b"")
        self.assertEqual(0, lookup[0]["offset"])
        self.assertEqual(fam0_count, lookup[0]["count"])

        # walk the lookup table and verify each group is family-homogeneous
        # and matches the offsets/counts partitioning
        pos = 0
        seen_fams = []
        for row in lookup:
            cnt = int(row["count"])
            group = reordered[pos : pos + cnt]
            if cnt > 0:
                fams_in_group = {(0 if hid == b"" else int(hid[4:11].decode())) for hid in group}
                self.assertEqual(1, len(fams_in_group), "group not homogeneous in family")
                seen_fams.append(fams_in_group.pop())
            pos += cnt
        self.assertEqual(len(hog_ids), pos)
        self.assertEqual(sorted(seen_fams), seen_fams)

    def test_rel_char_format(self):
        hog_ids = numpy.array(
            [
                b"HOG:A0000001.1a",
                b"HOG:A0000001.1b",
                b"HOG:A0000003",
            ],
            dtype="S255",
        )
        db = self._build_db_with_entries(hog_ids)
        db.add_protein_hog_ids(hog_ids)
        entries = db.h5.get_node("/Protein/Entries").read()
        numpy.testing.assert_array_equal(hog_ids, entries["OmaHOG"])

        lookup = db.h5.get_node("/Protein/FamIndex/Lookup").read()
        entry_idx = db.h5.get_node("/Protein/FamIndex/EntryIdx").read()
        reordered = hog_ids[entry_idx]
        pos = 0
        seen_fams = []
        for row in lookup:
            cnt = int(row["count"])
            group = reordered[pos : pos + cnt]
            if cnt > 0:
                # family 0 is always inserted as a leading placeholder row
                # (offset=0, count=0) when no entry actually has an empty
                # OmaHOG, so skip homogeneity checks on empty groups.
                fams_in_group = {int(hid[5:12].decode()) for hid in group}
                self.assertEqual(1, len(fams_in_group), "group not homogeneous in family")
                seen_fams.append(fams_in_group.pop())
            pos += cnt
        self.assertEqual(len(hog_ids), pos)
        self.assertEqual(sorted(seen_fams), seen_fams)


class IdentifyMainVariantsTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        path = os.path.join(self.tmpdir, "splice.h5")
        self.db = DBBuilder(path)
        self.db.__enter__()
        self.addCleanup(self.db.__exit__, None, None, None)
        self.vp_tab = self.db.h5.create_table(
            "/", "VPairs", description=tables.dtype_from_descr(tablefmt.PairwiseRelationTable)
        )

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _entries(self, n, **cols):
        e = numpy.zeros(n, dtype=tables.dtype_from_descr(tablefmt.ProteinTable))
        e["EntryNr"] = numpy.arange(1, n + 1)
        for k, v in cols.items():
            e[k] = v
        return e

    def _call(self, entries, offset=0, splice_size=None):
        splice_size = splice_size or (offset + len(entries))
        splice_arr = numpy.zeros(splice_size, dtype=numpy.int32)
        idx = list(range(1, len(entries) + 1))  # unused directly; real call builds idx from grp
        self.db._identify_main_variants(
            splice_groups=[list(range(1, len(entries) + 1))],
            splice_arr=splice_arr,
            entries=entries,
            offset=offset,
            vp_tab=self.vp_tab,
        )
        return splice_arr

    def test_single_oma_group_wins(self):
        # the "winning" EntryNr is broadcast to *all* members of the splice
        # group (this marks every splice variant's AltSpliceVariant with the
        # same "main isoform" EntryNr), not just written at the winner's own
        # position.
        entries = self._entries(3)
        entries["OmaGroup"][0] = 5
        splice_arr = self._call(entries)
        expected = entries["EntryNr"][0]
        numpy.testing.assert_array_equal([expected, expected, expected], splice_arr)

    def test_two_oma_group_raises(self):
        entries = self._entries(3)
        entries["OmaGroup"][0] = 5
        entries["OmaGroup"][1] = 6
        with self.assertRaises(DBConsistencyError):
            self._call(entries)

    def test_single_omahog_wins_when_no_omagroup(self):
        entries = self._entries(3)
        entries["OmaHOG"][1] = b"HOG:0000001"
        splice_arr = self._call(entries)
        expected = entries["EntryNr"][1]
        numpy.testing.assert_array_equal([expected, expected, expected], splice_arr)

    def test_two_omahog_raises(self):
        entries = self._entries(3)
        entries["OmaHOG"][0] = b"HOG:0000001"
        entries["OmaHOG"][1] = b"HOG:0000002"
        with self.assertRaises(DBConsistencyError):
            self._call(entries)

    def test_vp_choice_is_overwritten_by_longest_seq_BUG(self):
        # BUG-PIN: pyoma/browser/build/builder.py::DBBuilder._identify_main_variants
        # The "if len(vp) == 1: splice_arr[...] = ..." branch (numpy nonzero of
        # nr_vps) has no `continue`, unlike the OmaGroup/OmaHOG branches above
        # it. Execution therefore always falls through to the final
        # unconditional `splice_arr[idx + offset] = ent["EntryNr"][argmax(...)]`
        # line, silently discarding the pairwise-ortholog-based choice in favor
        # of "longest sequence wins" -- even when a unique variant with a
        # pairwise ortholog was found. This test pins the CURRENT (arguably
        # buggy) behavior; it does not validate it is correct.
        entries = self._entries(3, SeqBufferLength=[10, 999, 20])
        # give entry #1 (index 0) a pairwise ortholog
        self.vp_tab.append([(entries["EntryNr"][0], 42, 0, 0.0, 0.0, 0.0, 0.0, 0.0)])
        self.vp_tab.flush()
        splice_arr = self._call(entries)
        # If the vp-based choice were respected, every position would equal
        # entries["EntryNr"][0]. Instead, the longest SeqBufferLength
        # (index 1) wins everywhere.
        expected = entries["EntryNr"][1]
        numpy.testing.assert_array_equal([expected, expected, expected], splice_arr)

    def test_longest_seq_wins_with_no_group_hog_or_vp(self):
        entries = self._entries(3, SeqBufferLength=[10, 999, 20])
        splice_arr = self._call(entries)
        expected = entries["EntryNr"][1]
        numpy.testing.assert_array_equal([expected, expected, expected], splice_arr)

    def test_two_vp_raises(self):
        entries = self._entries(3, SeqBufferLength=[10, 20, 30])
        self.vp_tab.append(
            [
                (entries["EntryNr"][0], 42, 0, 0.0, 0.0, 0.0, 0.0, 0.0),
                (entries["EntryNr"][1], 43, 0, 0.0, 0.0, 0.0, 0.0, 0.0),
            ]
        )
        self.vp_tab.flush()
        with self.assertRaises(DBConsistencyError):
            self._call(entries)


class UpdateNrGenesTest(unittest.TestCase):
    def test_update_nr_genes(self):
        tmpdir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmpdir, ignore_errors=True)
        path = os.path.join(tmpdir, "nrgenes.h5")
        db = DBBuilder(path)
        db.__enter__()
        self.addCleanup(db.__exit__, None, None, None)

        genomes = numpy.zeros(2, dtype=tables.dtype_from_descr(tablefmt.GenomeTable))
        genomes[0]["UniProtSpeciesCode"] = b"HUMA"
        genomes[0]["EntryOff"] = 0
        genomes[0]["TotEntries"] = 4
        genomes[1]["UniProtSpeciesCode"] = b"MOUS"
        genomes[1]["EntryOff"] = 4
        genomes[1]["TotEntries"] = 3
        db.h5.create_table("/", "Genome", obj=genomes)

        entries = numpy.zeros(7, dtype=tables.dtype_from_descr(tablefmt.ProteinTable))
        entries["EntryNr"] = numpy.arange(1, 8)
        # genome 1 (entries 1-4): entry 2 is an alt splice variant of entry 1;
        # entries 3, 4 are main isoforms (AltSpliceVariant == 0) -> 3 genes
        entries["AltSpliceVariant"] = [1, 1, 0, 0, 0, 6, 6]
        # genome 2 (entries 5-7): entry 5 main (0), 6 main (==EntryNr), 7 alt splice of 6
        # -> 2 genes
        db.h5.create_group("/Protein", "Root", createparents=True)
        db.h5.create_table("/Protein", "Entries", obj=entries)

        db.update_nr_genes()

        genome_after = db.h5.get_node("/Genome").read()
        self.assertEqual(3, genome_after[0]["TotGenes"])
        self.assertEqual(2, genome_after[1]["TotGenes"])


class AddSequenceIndexTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _random_sequences(self, n, min_len=5, max_len=15, seed=123):
        rng = random.Random(seed)
        letters = [d.decode() for d in DIGITS_AA.tolist()]
        seqs = []
        for _ in range(n):
            length = rng.randint(min_len, max_len)
            seqs.append("".join(rng.choice(letters) for _ in range(length)))
        return seqs

    def _expected_kmer_entries(self, sequences, k):
        expected = {}
        for entry_nr, seq in enumerate(sequences, start=1):
            enc = KmerEncoder(k, is_protein=True)
            for code in enc.decompose(seq.encode("ascii")):
                expected.setdefault(code, []).append(entry_nr)
        return expected

    def _check_for_k(self, sequences, k):
        seqs = b"".join(s.encode("ascii") + b" " for s in sequences)
        path = os.path.join(self.tmpdir, f"sa_k{k}.h5")
        with DBBuilder(path) as db:
            db.add_sequence_index(seqs, nr_entries=len(sequences), k=k)

        with tables.open_file(path, mode="r") as h5:
            kmer_lookup = h5.root.Protein.KmerLookup
            expected = self._expected_kmer_entries(sequences, k)
            n_codes = KmerEncoder(k).n
            self.assertEqual(n_codes, len(kmer_lookup))

            total_found = 0
            for code in range(n_codes):
                got = sorted(int(x) for x in kmer_lookup[code])
                exp = sorted(expected.get(code, []))
                self.assertEqual(exp, got, f"mismatch at kmer code {code} for k={k}")
                total_found += len(got)

            total_expected_positions = sum(max(0, len(s) - k + 1) for s in sequences)
            self.assertEqual(total_expected_positions, total_found)

            keys = h5.root.Protein.SuffixArrayIndexKeys.read()
            pos = h5.root.Protein.SuffixArrayIndexPos.read()
            self.assertEqual(len(keys), len(pos))
            self.assertEqual(n_codes, keys[-1])
            self.assertEqual(len(seqs), pos[-1])

    # Documents why add_sequence_index (builder.py:816) needs `...[0][0]`
    # rather than `...[0]` to unpack `first_code`: kmer_codes_for_positions
    # returns a 2-tuple `(codes, good)`, and `codes` is itself a numpy array
    # (shape (1,) here, since `sa_f[0:1]` has length 1), not a scalar.
    # `int(<len-1 array>)` was only a (deprecated) implicit conversion under
    # numpy<2; under numpy>=2 it raises
    # `TypeError: only 0-dimensional arrays can be converted to Python
    # scalars`. Previously this made add_sequence_index fail unconditionally
    # whenever `L = len(sa_f) - k > 0` -- true for essentially any real,
    # non-trivial set of sequences. Fixed at builder.py:816 (`...[0][0]`).
    def test_kmer_codes_for_positions_returns_array_not_scalar(self):
        from pyoma.browser.build.builder import kmer_codes_for_positions

        seqs = b"ACDEFG ACDEFG "
        seqs_np = numpy.frombuffer(seqs, dtype=numpy.uint8)
        map256 = numpy.full(256, 255, dtype=numpy.uint8)
        for i, aa in enumerate(DIGITS_AA):
            map256[ord(aa)] = i
        k = 2
        dtype_sa = numpy.dtype("int64")
        sa_f = numpy.array([0, 1, 2], dtype=dtype_sa)
        codes_tuple_first = kmer_codes_for_positions(sa_f[0:1], k, seqs_np, dtype_sa, map256, len(DIGITS_AA))[0]
        self.assertEqual((1,), codes_tuple_first.shape)
        with self.assertRaises(TypeError):
            int(codes_tuple_first)

    def test_kmer_lookup_matches_brute_force_k2(self):
        sequences = self._random_sequences(15, seed=7)
        self._check_for_k(sequences, 2)

    def test_kmer_lookup_matches_brute_force_k3(self):
        sequences = self._random_sequences(15, seed=11)
        self._check_for_k(sequences, 3)

    def test_kmer_lookup_matches_brute_force_k4_small_buffer(self):
        # Regression test for a kmer-code integer overflow bug in
        # suffixarray_helper.kmer_codes_for_positions: it used to accumulate
        # kmer codes in `dtype_sa` (the suffix-array *position* dtype, sized
        # only for the sequence-buffer length -- PySAIS returns uint16 for any
        # buffer <~64KB). Kmer codes scale as alphabet_size**k independently of
        # buffer length: for k=4 with the 21-letter AA alphabet, codes go up to
        # 21**4 - 1 = 194480, which overflows uint16 (max 65535) and silently
        # wraps, corrupting the index -- but only when the sequence buffer is
        # small enough to get a narrow dtype_sa (real production genome
        # databases are always well above 64KB, so this never showed up there).
        # This test uses a small buffer (well under 64KB) with k=4 to pin the
        # fix (suffixarray_helper.dtype_for_kmer_codes, sized independently of
        # dtype_sa).
        sequences = self._random_sequences(30, seed=13)
        seqs = b"".join(s.encode("ascii") + b" " for s in sequences)
        self.assertLess(len(seqs), 65536, "test buffer must stay small to exercise a narrow dtype_sa")
        self._check_for_k(sequences, 4)


if __name__ == "__main__":
    unittest.main()
