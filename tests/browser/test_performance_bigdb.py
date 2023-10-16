from unittest import mock, SkipTest
import time
import os
from .test_db import TestWithDbInstance
from pyoma.browser.search import XRefSearch


class TestWithBigOmaDB(TestWithDbInstance):
    db_file_name = "All.Dec2021/data/OmaServer.h5"

    @classmethod
    def setUpClass(cls):
        with mock.patch.dict(
            os.environ, {"PYOMA_DB_PATH": "/Volumes/TOSHIBA EXT/browser"}
        ):
            try:
                super().setUpClass()
            except IOError:
                raise SkipTest(
                    "Big database '{}' is not available".format(cls.db_file_name)
                )


class XRefSearchTest(TestWithBigOmaDB):
    def test_existing_example(self):
        for query in ("TP53", "INS"):
            with self.subTest(query=query):
                s = XRefSearch(self.db, query, max_matches=80)
                hits = 0
                for p in s.search_entries():
                    for source, val in p.xref_data.items():
                        for el in val:
                            self.assertIn(query.lower(), el.lower())
                            hits += 1
                self.assertGreater(hits, 0)


if __name__ == "__main__":
    import unittest

    unittest.main()
