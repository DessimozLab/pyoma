from __future__ import unicode_literals
from future.standard_library import install_aliases

install_aliases()

import random
import re
import unittest

try:
    import unittest.mock as mock
except ImportError:
    import mock as mock
import numpy
import tables
import string
import itertools
from pyoma.browser import suffixsearch

CHARS = numpy.fromiter(
    itertools.chain(string.ascii_letters, string.digits, "-_."), dtype="S1"
)


class StringGenerator(object):
    def __init__(self):
        self.buffer = None
        self.pos = 0
        self.fill_buffer()

    def fill_buffer(self):
        self.buffer = numpy.random.choice(CHARS, 2 ** 20)
        self.pos = 0

    def get_string(self, length):
        if self.pos + length > len(self.buffer):
            self.fill_buffer()
        res = self.buffer[self.pos : self.pos + length].tostring()
        self.pos += length
        return res


def create_h5_with_table(nr_row):
    f = tables.open_file(
        "testsuffix.h5", "w", driver="H5FD_CORE", driver_core_backing_store=0
    )
    tab_def = numpy.dtype(
        [("CharCol50", "S50"), ("DblCol", "f8"), ("CharCol200", "S200")]
    )
    tab = numpy.zeros(nr_row, tab_def)
    generator = StringGenerator()
    tab["CharCol50"] = [
        generator.get_string(l) for l in numpy.random.randint(0, 51, size=nr_row)
    ]
    tab["DblCol"] = numpy.random.rand(nr_row)
    tab["CharCol200"] = [
        generator.get_string(l) for l in numpy.random.randint(0, 201, size=nr_row)
    ]
    f.create_table("/test", "table", obj=tab, createparents=True)
    return f


def create_h5_with_varlen_string_col(nr_row, off_col_type):
    f = tables.open_file(
        "testsuffix.h5", "w", driver="H5FD_CORE", driver_core_backing_store=0
    )
    tab_def = numpy.dtype(
        [("CharCol", "S800"), ("VarCharOff", off_col_type), ("VarCharLen", "i4")]
    )
    buf = f.create_earray(
        "/test", "buffer", tables.StringAtom(1), (0,), createparents=True
    )
    tab = numpy.zeros(nr_row + 1, tab_def)
    generator = StringGenerator()
    off = 0
    for i in range(nr_row):
        l = random.randint(0, 800)
        val = numpy.ndarray(
            (l,), buffer=generator.get_string(l), dtype=tables.StringAtom(1)
        )
        buf.append(val)
        tab["VarCharOff"][i] = off
        tab["VarCharLen"][i] = l
        tab["CharCol"][i] = val.tostring()
        off += l
    # last row does not contain any data, simulate not set values in Offset array
    f.create_table("/test", "table", obj=tab)
    return f


class SuffixArraySearchTests(unittest.TestCase):
    def setUp(self):
        self.h5 = create_h5_with_table(5000)

    def tearDown(self):
        self.h5.close()

    def test_create_index_for_CharCol50(self):
        suf = suffixsearch.SuffixIndexBuilderStringCol(
            self.h5.get_node("/test/table"), "CharCol50", "/test", ignore_case=True
        )
        suf()
        self.assertEqual(
            self.h5.get_node_attr("/test/table", "CharCol50_suffixindexnode"), "/test"
        )
        self.assertEqual(
            self.h5.get_node_attr("/test", "CharCol50_buffer"), "/test/CharCol50_buffer"
        )
        self.assertEqual(
            self.h5.get_node_attr("/test", "CharCol50_offset"), "/test/CharCol50_offset"
        )
        self.assertEqual(self.h5.get_node_attr("/test", "CharCol50_ignore_case"), True)
        try:
            suf = self.h5.get_node("/test/CharCol50_suffix")
            buf = self.h5.get_node("/test/CharCol50_buffer")
            off = self.h5.get_node("/test/CharCol50_offset")
        except tables.NoSuchNodeError as e:
            self.assertTrue(False, "Suffix index data not complete: {}".format(e))

    def test_create_index_more_than_one_column(self):
        cols = ("CharCol50", "CharCol200")
        for col in cols:
            suf = suffixsearch.SuffixIndexBuilderStringCol(
                self.h5.get_node("/test/table"), col, "/test", ignore_case=True
            )
            suf()
        for col in cols:
            self.assertEqual(
                "/test/{}_buffer".format(col),
                self.h5.get_node_attr("/test", col + "_buffer"),
            )
            self.assertEqual(
                "/test/{}_offset".format(col),
                self.h5.get_node_attr("/test", col + "_offset"),
            )

    def test_search_case_insensitive(self):
        tab = self.h5.get_node("/test/table")
        suffixsearch.SuffixIndexBuilderStringCol(
            tab, "CharCol50", "/test", ignore_case=True
        )()
        search = suffixsearch.SuffixSearcher.from_tablecolumn(tab, "CharCol50")
        target_row = numpy.random.randint(0, len(tab), size=200)
        for k, full_query in enumerate(
            tab.read_coordinates(target_row, field="CharCol50")
        ):
            if len(full_query) < 2:
                continue
            qlen = random.randint(2, len(full_query) + 1)
            strt = random.randint(0, len(full_query) + 1 - qlen)
            query = full_query[strt : strt + qlen]
            res = search.find(query)
            self.assertIn(
                target_row[k],
                res,
                "could not find {} (full value {}, row {}) in results ({})".format(
                    query, full_query, target_row[k], res
                ),
            )
            for target in res:
                self.assertIsNotNone(
                    re.search(query, tab[target]["CharCol50"], re.IGNORECASE)
                )

    def test_search_string_instance(self):
        tab = self.h5.get_node("/test/table")
        suffixsearch.SuffixIndexBuilderStringCol(
            tab, "CharCol50", "/test", ignore_case=True
        )()
        search = suffixsearch.SuffixSearcher.from_tablecolumn(tab, "CharCol50")
        for trial in range(100):
            row_nr = random.randint(0, len(tab))
            query = tab[row_nr]["CharCol50"].decode()
            if len(query) > 3:
                break
        else:
            self.assertTrue(False, "table has no valid row to search")
        res = search.find(query[2:])
        self.assertIn(row_nr, res)

    def test_raises_on_non_char_column(self):
        tab = self.h5.get_node("/test/table")
        with self.assertRaises(TypeError):
            suffixsearch.SuffixIndexBuilderStringCol(
                tab, "DblCol", "/test", ignore_case=True
            )()


class SuffixArrayVarLenSearchTestsCaseInsensitive(unittest.TestCase):
    ignore_case = True

    def setUp(self):
        self.h5 = create_h5_with_varlen_string_col(5000, numpy.int32)
        self.build_index()

    def build_index(self):
        tab = self.h5.get_node("/test/table")
        suffixsearch.SuffixIndexBuilderVarStringCol(
            tab,
            "VarCharOff",
            self.h5.get_node("/test/buffer"),
            "/test",
            ignore_case=self.ignore_case,
        )()

    def tearDown(self):
        self.h5.close()

    def test_var_and_fixed_col_have_same_values(self):
        tab = self.h5.get_node("/test/table")
        buf = self.h5.get_node("/test/buffer")
        for row in tab:
            from_buf = buf[
                row["VarCharOff"] : (row["VarCharOff"] + row["VarCharLen"])
            ].tostring()
            self.assertEqual(row["CharCol"], from_buf)

    def test_existance_of_auxillary_buffers(self):
        self.assertEqual(
            self.h5.get_node_attr("/test/table", "VarCharOff_suffixindexnode"), "/test"
        )
        if self.ignore_case:
            self.assertEqual(
                self.h5.get_node_attr("/test", "VarCharOff_buffer"),
                "/test/VarCharOff_buffer",
            )
        self.assertEqual(
            self.h5.get_node_attr("/test", "VarCharOff_offset"),
            "/test/VarCharOff_offset",
        )
        self.assertEqual(
            self.h5.get_node_attr("/test", "VarCharOff_ignore_case"), self.ignore_case
        )
        try:
            suf = self.h5.get_node("/test/VarCharOff_suffix")
            off = self.h5.get_node("/test/VarCharOff_offset")
            if self.ignore_case:
                buf = self.h5.get_node("/test/VarCharOff_buffer")
            else:
                with self.assertRaises(tables.NoSuchNodeError):
                    self.h5.get_node("/test/VarCharOff_buffer")

        except tables.NoSuchNodeError as e:
            self.assertTrue(False, "Suffix index data not complete: {}".format(e))

    def test_search_pattern(self):
        tab = self.h5.get_node("/test/table")
        search = suffixsearch.SuffixSearcher.from_tablecolumn(tab, "VarCharOff")
        target_row = numpy.random.randint(0, len(tab), size=200)
        for k, full_query in enumerate(
            tab.read_coordinates(target_row, field="CharCol")
        ):
            if len(full_query) < 2:
                continue
            qlen = random.randint(2, len(full_query) + 1)
            strt = random.randint(0, len(full_query) + 1 - qlen)
            query = full_query[strt : strt + qlen]
            res = search.find(query)
            self.assertIn(
                target_row[k],
                res,
                "could not find {} (full value {}, row {}) in results ({})".format(
                    query, full_query, target_row[k], res
                ),
            )
            for target in res:
                if self.ignore_case:
                    self.assertIsNotNone(
                        re.search(query, tab[target]["CharCol"], flags=re.IGNORECASE),
                        "Query {} not found in {}. target: {}".format(
                            query, tab[target]["CharCol"], target
                        ),
                    )
                else:
                    if query not in tab[target]["CharCol"]:
                        # we didn't add a sentinel, could be a combination with the next value
                        combined = b"".join(tab[target : target + 2]["CharCol"])
                        self.assertIn(query, combined)

    def test_non_existing_pattern(self):
        tab = self.h5.get_node("/test/table")
        search = suffixsearch.SuffixSearcher.from_tablecolumn(tab, "VarCharOff")
        res = search.find("ZZZZzz")
        self.assertEqual(0, len(res))


class SuffixArrayVarLenSearchTestsCaseSensitive(
    SuffixArrayVarLenSearchTestsCaseInsensitive
):
    ignore_case = False


class SuffixBuilderFactoryTester(unittest.TestCase):
    offset_col_type = numpy.int16

    def setUp(self):
        self.h5 = create_h5_with_varlen_string_col(100, self.offset_col_type)

    def tearDown(self):
        self.h5.close()

    @mock.patch("pyoma.browser.suffixsearch.SuffixIndexBuilderVarStringCol")
    @mock.patch("pyoma.browser.suffixsearch.SuffixIndexBuilderStringCol")
    def test_fix_string_column(self, FixStrSuffixMock, VarStrSuffixMock):
        tab = self.h5.get_node("/test/table")
        suffixsearch.create_suffix_index(tab, "CharCol", ignore_case=True)
        FixStrSuffixMock.assert_called_once_with(
            tab, "CharCol", self.h5.get_node("/test/_si_table"), ignore_case=True
        )
        VarStrSuffixMock.assert_not_called()

    @mock.patch("pyoma.browser.suffixsearch.SuffixIndexBuilderVarStringCol")
    @mock.patch("pyoma.browser.suffixsearch.SuffixIndexBuilderStringCol")
    def test_var_string_column(self, FixStrSuffixMock, VarStrSuffixMock):
        tab = self.h5.get_node("/test/table")
        buf = "/test/buffer"
        suffixsearch.create_suffix_index(tab, "VarCharOff", buf, ignore_case=True)
        FixStrSuffixMock.assert_not_called()
        VarStrSuffixMock.called_once_with(
            tab,
            "VarCharOff",
            self.h5.get_node(buf),
            self.h5.get_node("/test/_si_table"),
            ignore_case=True,
        )


class SuffixBuilderFactoryTesterUnsignedColOffset(SuffixBuilderFactoryTester):
    offset_col_type = numpy.uint64
