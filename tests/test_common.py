import os
import tempfile
import unittest
import tempfile
import gzip
import bz2
from pyoma.common import auto_open


class AutoOpenBaseTest(unittest.TestCase):
    filesuffix = None

    def setUp(self) -> None:
        self.expected_text = """This is a test text. Let's see if we can properly load it."""
        with tempfile.NamedTemporaryFile(suffix=self.filesuffix, delete=False) as fh:
            self.testfilename = fh.name

    def tearDown(self) -> None:
        os.remove(self.testfilename)


class AutoOpenRegularReadTest(AutoOpenBaseTest):
    def setUp(self) -> None:
        super().setUp()
        self.store_text(self.testfilename)

    def store_text(self, fn):
        with open(fn, mode="wt", encoding="utf-8") as fh:
            fh.write(self.expected_text)

    def test_read_as_text(self):
        with auto_open(self.testfilename, "rt") as fh:
            res = fh.read()
        self.assertEqual(res, self.expected_text)

    def test_read_as_bytes(self):
        with auto_open(self.testfilename, "rb") as fh:
            res = fh.read()
        self.assertEqual(res, self.expected_text.encode("utf-8"))


class AutoOpenGzipTest(AutoOpenRegularReadTest):
    def store_text(self, fn):
        with gzip.open(fn, mode="wt", encoding="utf-8") as fh:
            fh.write(self.expected_text)


class AutoOpenBz2Test(AutoOpenRegularReadTest):
    def store_text(self, fn):
        with bz2.open(fn, mode="wt", encoding="utf-8") as fh:
            fh.write(self.expected_text)


class AutoOpenRegularWriteTest(AutoOpenBaseTest):
    def expected_text_start(self):
        return self.expected_text.encode("utf-8")

    def test_write(self):
        with auto_open(self.testfilename, "wt") as fh:
            fh.write(self.expected_text)
        with open(self.testfilename, "rb") as fh:
            res = fh.read()
            self.assertTrue(res.startswith(self.expected_text_start()))


class AutoOpenGzipWriteTest(AutoOpenRegularWriteTest):
    filesuffix = ".gz"

    def expected_text_start(self):
        return b"\x1f\x8b\x08"


class AutoOpenBz2WriteTest(AutoOpenRegularWriteTest):
    filesuffix = ".bz2"

    def expected_text_start(self):
        return b"\x42\x5a\x68"


if __name__ == "__main__":
    unittest.main()
