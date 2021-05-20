import logging
import gzip
import bz2
import os
from io import BytesIO, StringIO

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

    # File opening. This is based on the example on SO here:
    # http://stackoverflow.com/a/26986344
    fmagic = {b"\x1f\x8b\x08": gzip.open, b"\x42\x5a\x68": bz2.open}

    if os.path.isfile(fn) and os.stat(fn).st_size > 0:
        with open(fn, "rb") as fp:
            fs = fp.read(max([len(x) for x in fmagic]))
        for (magic, _open) in fmagic.items():
            if fs.startswith(magic):
                return _open(fn, *args, **kwargs)
    else:
        if fn.endswith("gz"):
            return gzip.open(fn, *args, **kwargs)
        elif fn.endswith("bz2"):
            return bz2.open(fn, *args, **kwargs)
    return open(fn, *args, **kwargs)
