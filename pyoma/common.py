import collections
import itertools
import logging
import gzip
import bz2
import os
from io import BytesIO, StringIO
from pathlib import Path
import lzma

package_logger = logging.getLogger("pyoma")
package_logger.addHandler(logging.NullHandler())


def auto_open(fn, *args, **kwargs):
    """function to open regular or compressed files for read / write.

    This function opens files based on their "magic bytes". Supports bz2
    and gzip. If it finds neither of these, presumption is it is a
    standard, uncompressed file.

    Example::

        with auto_open("/path/to/file/maybe/compressed", mode="rb") as fh:
            fh.read()

        with auto_open("/tmp/test.txt.gz", mode="wb") as fh:
            fh.write("my big testfile")

    :param fn: either a string of an existing or new file path, or
        a BytesIO handle
    :param **kwargs: additional arguments that are understood by the
        underlying open handler
    :returns: a file handler
    """
    if isinstance(fn, (BytesIO, StringIO)):
        return fn
    if isinstance(fn, (Path, os.DirEntry)):
        fn = str(fn)

    # File opening. This is based on the example on SO here:
    # http://stackoverflow.com/a/26986344
    fmagic = {
        b"\x1f\x8b\x08": gzip.open,
        b"\x42\x5a\x68": bz2.open,
        b"\xfd\x37\x7a\x58\x5a\x00": lzma.open,  # xz
    }

    # Try detecting via magic bytes (if file exists and non-empty)
    if os.path.isfile(fn) and os.stat(fn).st_size > 0:
        with open(fn, "rb") as fp:
            fs = fp.read(max(map(len, fmagic)))
        for magic, _open in fmagic.items():
            if fs.startswith(magic):
                return _open(fn, *args, **kwargs)

    # Fallback to detection via file extension
    suffix = os.path.splitext(fn)[-1].lower()
    if suffix == ".gz":
        return gzip.open(fn, *args, **kwargs)
    elif suffix == ".bz2":
        return bz2.open(fn, *args, **kwargs)
    elif suffix == ".xz":
        return lzma.open(fn, *args, **kwargs)
    return open(fn, *args, **kwargs)


def count_elements(iterable):
    """return the number of elements in an iterator in the most efficient way.

    Be aware that for unbound iterators, this method won't terminate!
    :param iterable: an iterable object.
    """
    counter = itertools.count()
    collections.deque(zip(iterable, counter), maxlen=0)  # (consume at C speed)
    return next(counter)
