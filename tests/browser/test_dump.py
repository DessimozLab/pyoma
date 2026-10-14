import io
import os
import tempfile
import unittest

import numpy
import tables

from pyoma.browser.build.dump import _write_roothog_descriptions
from pyoma.browser.tablefmt import RootHOGMetaTable


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
