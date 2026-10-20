import io
import os
import tempfile
import unittest

import numpy
import tables

from pyoma.browser.build.dump import _write_roothog_descriptions, _write_xrefs
from pyoma.browser.tablefmt import RootHOGMetaTable, XRefTable


class RootHOGDescriptionsDumpTest(unittest.TestCase):
    def setUp(self):
        fd, self.path = tempfile.mkstemp(suffix=".h5")
        os.close(fd)
        keywords = {1: "Kinase, ATP binding", 2: "", 3: "Récepteur"}
        buf, rows = b"", []
        for fam in (3, 1, 2):  # deliberately unordered
            kw = keywords[fam].encode("utf-8")
            rows.append((fam, len(buf), len(kw), 0, 0))
            buf += kw
        with tables.open_file(self.path, "w") as h5:
            tab = h5.create_table("/RootHOG", "MetaData", RootHOGMetaTable, createparents=True)
            tab.append(rows)
            kwbuf = h5.create_earray("/RootHOG", "KeywordBuffer", tables.StringAtom(1), (0,))
            kwbuf.append(numpy.frombuffer(buf, dtype="S1"))

    def tearDown(self):
        os.remove(self.path)

    def test_writes_sorted_tsv(self):
        fh = io.StringIO()
        with tables.open_file(self.path) as h5:
            n = _write_roothog_descriptions(fh, h5, release_char="A")
        self.assertEqual(3, n)
        self.assertEqual(
            "HOG:A0000001\tKinase, ATP binding\nHOG:A0000002\t\nHOG:A0000003\tRécepteur\n",
            fh.getvalue(),
        )


if __name__ == "__main__":
    unittest.main()


class XRefsDumpTest(unittest.TestCase):
    def setUp(self):
        fd, self.path = tempfile.mkstemp(suffix=".h5")
        os.close(fd)
        self.rows = [
            (1, "SourceAC", "AC1"),
            (1, "UniProtKB/SwissProt", "P12345"),
            (1, "SourceID", "ID1"),
            (1, "SourceAC", "AC1-2"),
            (1, "SourceAC", "AC1"),
            (2, "UniProtKB/SwissProt", "P99999"),
            (3, "SourceID", "ID3"),
            (4, "SourceAC", "AC4"),
        ]
        self.h5 = tables.open_file(self.path, "w")
        self.tab = self.h5.create_table("/", "XRef", XRefTable)
        self.set_rows(self.rows)

    def set_rows(self, rows):
        enum = self.tab.get_enum("XRefSource")
        self.tab.remove_rows(0, self.tab.nrows)
        for enr, src, xref in rows:
            self.tab.row["EntryNr"] = enr
            self.tab.row["XRefSource"] = enum[src]
            self.tab.row["XRefId"] = xref.encode("utf-8")
            self.tab.row.append()
        self.tab.flush()

    def tearDown(self):
        self.h5.close()
        os.remove(self.path)

    def dump(self, sources):
        handles = {}
        outputs = {src: handles.setdefault(name, io.StringIO()) for src, name in sources.items()}
        _write_xrefs(self.tab, outputs, lambda enr: f"OMA{enr:05d}")
        return {name: fh.getvalue() for name, fh in handles.items()}

    def test_source_id_and_ac_on_one_line(self):
        res = self.dump({"SourceID": "source", "SourceAC": "source", "UniProtKB/SwissProt": "uniprot"})
        self.assertEqual("OMA00001\tID1\tAC1;AC1-2\nOMA00003\tID3\t\nOMA00004\t\tAC4\n", res["source"])
        self.assertEqual("OMA00001\tP12345\nOMA00002\tP99999\n", res["uniprot"])

    def test_unsorted_table_raises_when_combining(self):
        self.set_rows(self.rows[5:] + self.rows[:5])
        with self.assertRaises(ValueError):
            self.dump({"SourceID": "source", "SourceAC": "source"})

    def test_unsorted_table_ok_without_source_output(self):
        self.set_rows(self.rows[5:] + self.rows[:5])
        res = self.dump({"UniProtKB/SwissProt": "uniprot"})
        self.assertEqual("OMA00002\tP99999\nOMA00001\tP12345\n", res["uniprot"])
