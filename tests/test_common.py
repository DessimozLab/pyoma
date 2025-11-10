import io
import lzma
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


class AutoOpenXzTest(AutoOpenRegularReadTest):
    filesuffix = ".xz"

    def store_text(self, fn):
        with lzma.open(fn, mode="wt", encoding="utf-8") as fh:
            fh.write(self.expected_text)


class AutoOpenBytesIOTest(unittest.TestCase):
    def test_bytesio_returns_itself(self):
        data = b"Hello World"
        buf = io.BytesIO(data)
        fh = auto_open(buf)
        self.assertIs(fh, buf)
        self.assertEqual(fh.read(), data)


class AutoOpenPathTest(AutoOpenRegularReadTest):
    def store_text(self, fn):
        with open(fn, "wt", encoding="utf-8") as fh:
            fh.write(self.expected_text)

    def test_pathlib_path(self):
        from pathlib import Path

        path_obj = Path(self.testfilename)
        with auto_open(path_obj, "rt") as fh:
            self.assertEqual(fh.read(), self.expected_text)


if __name__ == "__main__":
    unittest.main()
