from __future__ import division, print_function, absolute_import, unicode_literals

import time
from builtins import bytes, str, range, str
from bisect import bisect_left
import os
import numpy
import tables
from .models import KeyWrapper
import logging

logger = logging.getLogger(__name__)


class SuffixIndexBuilderStringCol(object):
    MEM_PER_CHUNK = 50 * 2 ** 20

    def __init__(self, tab, col, index_group, ignore_case):
        if not isinstance(tab, tables.Table):
            raise ValueError("tab argument must be a tables.Table object")
        if col not in tab.colnames:
            raise ValueError("column '{}' does not exist in {}".format(col, tab))
        self.h5 = tab._v_file
        self.tab = tab
        self.col = col
        self.index_group = (
            index_group
            if isinstance(index_group, tables.Group)
            else self.h5.get_node(index_group)
        )
        self.ignore_case = ignore_case

    def check_column_types_or_raise(self):
        if not numpy.issubdtype(self.tab.coldtypes[self.col].type, numpy.bytes_):
            raise TypeError(
                "column '{}' must be a character type column".format(self.col)
            )

    def get_expected_index_length(self):
        return self.tab.coldtypes[self.col].itemsize * len(self.tab)

    def _arrayname(self, kind):
        return self.col + "_" + kind

    def _remove_aux_arrays(self):
        for kind in ("suffix", "buffer", "offset"):
            try:
                n = self.get_aux_array_handle(kind)
                index_group_path = self.index_group._v_pathname
                if (
                    os.path.commonprefix([index_group_path, n._v_pathname])
                    == index_group_path
                ):
                    # auxilary buffer hangs on the index_group. we can remove it
                    logger.info(
                        "removing existing auxillary node: {}".format(n._v_pathname)
                    )
                    n.remove()
                    try:
                        self.h5.del_node_attr(self.index_group, self._arrayname(kind))
                    except AttributeError:
                        pass
            except tables.NoSuchNodeError:
                pass

    def create_aux_arrays(self):
        for kind, typ in (
            ("buffer", tables.StringAtom(1)),
            ("offset", tables.UInt32Atom()),
        ):
            arr_name = self._arrayname(kind)
            exp_rows = (
                self.get_expected_index_length() if kind == "buffer" else len(self.tab)
            )
            arr = self.h5.create_earray(
                self.index_group, arr_name, typ, (0,), expectedrows=exp_rows
            )
            self.h5.set_node_attr(self.index_group, arr_name, arr._v_pathname)

    def get_aux_array_handle(self, kind):
        if kind not in ("buffer", "offset", "suffix"):
            raise ValueError("Not a valid handle for suffix index")
        attr = self._arrayname(kind)
        try:
            path = self.h5.get_node_attr(self.index_group, attr)
        except AttributeError:
            raise tables.NoSuchNodeError("'{}' does not exist".format(attr))
        return self.h5.get_node(path)

    def build_index_buffer(self):
        chunksize = int(self.MEM_PER_CHUNK / self.tab.coldtypes[self.col].itemsize)
        total_offset = 0
        buffer, offset = (
            self.get_aux_array_handle(kind) for kind in ("buffer", "offset")
        )
        for chunk_start in range(0, len(self.tab), chunksize):
            # load fixed width string col in numpy array, compute actual string lengths and build
            # a \x00 delimited buffer of all values
            col_data = self.tab.read(
                start=chunk_start, stop=chunk_start + chunksize, field=self.col
            )
            data_lens = numpy.char.str_len(col_data)
            data_as_long_arr = col_data.view("S1")
            tot_len = (data_lens + 1).sum()
            buf = numpy.zeros(tot_len, dtype="S1")
            t = (data_lens + 1).cumsum()
            starts = numpy.roll(t, 1)
            starts[0] = 0
            ends = t - 1
            stride = col_data.strides[0]
            for i in range(len(col_data)):
                buf[starts[i] : ends[i]] = data_as_long_arr[
                    (stride * i) : (stride * i + data_lens[i])
                ]
            if self.ignore_case:
                buf = numpy.char.lower(buf)
            # update global offset and append chunked buffer
            starts += total_offset
            buffer.append(buf)
            offset.append(starts)
            total_offset += tot_len

    def build_suffix_array(self):
        from PySAIS import sais

        data = self.get_aux_array_handle("buffer")[:]
        suffix = sais(data)
        arr = self.h5.create_carray(
            self.index_group, self._arrayname("suffix"), obj=suffix
        )
        self.h5.set_node_attr(
            self.index_group, self._arrayname("suffix"), arr._v_pathname
        )
        self.h5.set_node_attr(
            self.index_group, self.col + "_ignore_case", self.ignore_case
        )
        self.h5.set_node_attr(
            self.tab, self.col + "_suffixindexnode", self.index_group._v_pathname
        )

    def __call__(self, force=True):
        self.check_column_types_or_raise()
        if force:
            self._remove_aux_arrays()
        self.create_aux_arrays()
        self.build_index_buffer()
        self.build_suffix_array()


class SuffixIndexBuilderVarStringCol(SuffixIndexBuilderStringCol):
    def __init__(self, tab, col, buffer, index_group, ignore_case):
        if not isinstance(buffer, tables.CArray):
            raise TypeError("buffer argument must be a tables.CArray instance")
        super(SuffixIndexBuilderVarStringCol, self).__init__(
            tab, col, index_group, ignore_case
        )
        self.orig_buffer = buffer

    def check_column_types_or_raise(self):
        if not numpy.issubdtype(self.orig_buffer.dtype.type, numpy.bytes_):
            raise TypeError(
                "buffer '{}' must be a character type column".format(
                    self.orig_buffer._v_pathname
                )
            )
        if not numpy.issubdtype(self.tab.coldtypes[self.col].type, numpy.integer):
            raise TypeError(
                "column '{}' must be an integer argument with offsets into the '{}' buffer array".format(
                    self.col, self.orig_buffer._v_pathname
                )
            )

    def get_expected_index_length(self):
        return len(self.orig_buffer) + len(
            self.tab
        )  # extra char per row for delimiting

    def create_aux_arrays(self):
        for kind, typ in (
            ("buffer", tables.StringAtom(1)),
            ("offset", tables.UInt32Atom()),
        ):
            arr_name = self._arrayname(kind)
            if not self.ignore_case and kind == "buffer":
                self.h5.set_node_attr(
                    self.index_group, arr_name, self.orig_buffer._v_pathname
                )
            else:
                exp_rows = (
                    self.get_expected_index_length()
                    if kind == "buffer"
                    else len(self.tab)
                )
                arr = self.h5.create_earray(
                    self.index_group, arr_name, typ, (0,), expectedrows=exp_rows
                )
                self.h5.set_node_attr(self.index_group, arr_name, arr._v_pathname)

    def build_index_buffer(self):
        off_data = self.h5.get_node(self.tab).read(field=self.col)
        # first, we need to fix cases where off_data is 0 in the middle of an array,
        # simply because it has not been set as its length is also 0. We do this
        # by setting the offsets to the next starting offset of a non-zero length
        # entry (by copying the next value from the back)
        # set every zero position except [0] to the one following in the buffer,
        # in case the buffer ends with an uninitialized value, set it to the
        # length of the entire buffer.
        if off_data[-1] == 0:
            off_data[-1] = len(self.orig_buffer)
        zero_positions = numpy.where(off_data == 0)[0]
        for idx in zero_positions[:0:-1]:
            off_data[idx] = off_data[idx + 1]

        if self.ignore_case:
            starts = off_data
            # now, set the stop position by using the next start postion (shift array by 1)
            stops = numpy.roll(off_data, -1)
            stops[-1] = len(self.orig_buffer)
            buf_arr = self.get_aux_array_handle("buffer")
            CHUNKSIZE = 2 ** 22  # 4MB
            tmp_buf = numpy.zeros(CHUNKSIZE, dtype="S1")
            tmp_buf_idx = 0
            for i in range(len(off_data)):
                if stops[i] - starts[i] + tmp_buf_idx >= CHUNKSIZE:
                    buf_arr.append(numpy.char.lower(tmp_buf[0:tmp_buf_idx]))
                    tmp_buf = numpy.zeros(CHUNKSIZE, dtype="S1")
                    tmp_buf_idx = 0
                tmp_buf[
                    tmp_buf_idx : (tmp_buf_idx + stops[i] - starts[i])
                ] = self.orig_buffer[starts[i] : stops[i]]
                tmp_buf_idx += (
                    stops[i] - starts[i] + 1
                )  # keep one extra '\x00' as separator
            if tmp_buf_idx > 0:
                buf_arr.append(numpy.char.lower(tmp_buf[0:tmp_buf_idx]))
            # add one extra position to each string that contains a '\x00'
            off_data += numpy.arange(0, len(off_data), 1, dtype=off_data.dtype)

        off_arr = self.get_aux_array_handle("offset")
        off_arr.append(off_data)


def create_or_load_index_group(tab, index_group=None):
    f = tab._v_file
    if index_group is None:
        idx_prefix, idx_name = tab._v_parent, "_si_" + tab._v_name
    else:
        idx_prefix, idx_name = index_group.rsplit("/", 1)

    try:
        return f.get_node(idx_prefix, name=idx_name)
    except tables.NoSuchNodeError:
        return f.create_group(idx_prefix, idx_name)


def create_suffix_index(
    tab, col, buffer=None, index_group=None, ignore_case=True, force=False
):
    """Create a suffix array from the given table column (or buffer).

    This function creates for a given fix-width char column a suffix array for
    fast lookup of (partial) matches. In this case, the `buffer` argument
    should be None. Alternatively, if a index should be build from a variable
    length column, i.e. a buffer array with offset and length attributes in
    the table, one should pass the offset column as `col` argument and the
    char buffer as `buffer` argument.

    By default, the class will create an array named <colname>_suffix,
    <colname>_buffer and <colname>_offset in a group node that is named and
    located either according to the argument `index_group` or as a sibling
    to `tab` named <tab>_Indexes.

    `ignore_case` will indicate that the index should be build case insensitive.
    In case an index is build from a variable length column (see above what
    we understand this is in this context) and `ignore_case` is False, the
    buffer and offset array will not be built again.

    :param tab: a handle to a h5 table node
    :param str col: the name the column to be indexed
    :param buffer: handle to buffer array
    :param str index_group: path to existing or non-existing group node
        where index data will be stored.
    :param bool ignore_case: case insensitive lookup
    :param force: overwrite existing suffix array data"""
    if not isinstance(tab, tables.Table):
        raise TypeError(
            "tab arguments must be a tables.Table (is a {})".format(type(tab))
        )
    if not isinstance(col, str) or col not in tab.colnames:
        raise ValueError("table {} does not have a column named '{}'".format(tab, col))
    typ = tab.coldtypes[col]
    if numpy.issubdtype(typ.type, numpy.bytes_):
        if typ.itemsize < 3:
            raise TypeError(
                "suffix arrays should be created for longer string columns only"
            )
        suffixbuilder = SuffixIndexBuilderStringCol(
            tab,
            col,
            create_or_load_index_group(tab, index_group),
            ignore_case=ignore_case,
        )
    elif numpy.issubdtype(typ.type, numpy.integer):
        if buffer is None:
            raise ValueError("buffer array must be specified for numeric table columns")
        if isinstance(buffer, str):
            buffer = tab._v_file.get_node(buffer)
        if not (
            isinstance(buffer, tables.CArray)
            and numpy.issubdtype(buffer, numpy.bytes_)
            and buffer.dtype.itemsize == 1
        ):
            raise TypeError(
                "buffer array must by a character tables.CArray instance (itemsize == 1)"
            )
        suffixbuilder = SuffixIndexBuilderVarStringCol(
            tab,
            col,
            buffer,
            create_or_load_index_group(tab, index_group),
            ignore_case=ignore_case,
        )
    else:
        raise TypeError("unsupported type of column to index")
    return suffixbuilder()


class SuffixSearcher(object):
    @classmethod
    def from_tablecolumn(cls, table, column, ignore_case=False):
        h5 = table._v_file
        try:
            idx_node = h5.get_node_attr(table, column + "_suffixindexnode")
            idx_node = h5.get_node(idx_node)
        except (AttributeError, tables.NoSuchNodeError) as e:
            raise SuffixIndexError(
                "Column {} of table {} does not seem to have a suffix index".format(
                    column, table
                )
            )
        try:
            suffix_arr = h5.get_node(idx_node, column + "_suffix")
            buffer_arr = h5.get_node(h5.get_node_attr(idx_node, column + "_buffer"))
            offset_arr = h5.get_node(h5.get_node_attr(idx_node, column + "_offset"))
        except (tables.NoSuchNodeError, AttributeError) as e:
            raise SuffixIndexInconsitency(
                "not all suffix index elements available: {}".format(e)
            )
        try:
            ignore_case = bool(h5.get_node_attr(idx_node, column + "_ignore_case"))
        except AttributeError:
            ignore_case = False
        return cls(suffix_arr, buffer_arr, offset_arr, ignore_case)

    @classmethod
    def from_index_node(cls, index_node, buffer=None, offset=None):
        import warnings

        warnings.warn(
            "initializing SuffixSearcher this way is deprecated. Use from_tablecolumn(table, columnname) instead"
        )
        if not isinstance(index_node, tables.Group):
            raise TypeError("expected a tables.Group node pointing to the index")
        try:
            buffer_arr = buffer if buffer else index_node._f_get_child("buffer")
            suffix_arr = index_node._f_get_child("suffix")
            offset_arr = offset if offset else index_node._f_get_child("offset")
            return cls(suffix_arr, buffer_arr, offset_arr, ignore_case=True)
        except tables.NoSuchNodeError as e:
            raise SuffixIndexInconsitency("suffix array data missing: {}".formate(e))

    def __init__(self, suffix, buffer, offset, ignore_case):
        self.suffix_arr = suffix
        self.buffer_arr = buffer
        self.offset_arr = offset[:]
        self.ignore_case = bool(ignore_case)

    def find(self, query):
        if isinstance(query, str):
            query = query.encode("utf-8")
        if self.ignore_case:
            query = query.lower()
        n = len(query)
        if n > 0:
            slicer = KeyWrapper(
                self.suffix_arr, key=lambda i: self.buffer_arr[i : (i + n)].tobytes()
            )
            t0 = time.time()
            ii = bisect_left(slicer, query)
            t1 = time.time()
            if ii and ii < len(self.suffix_arr) and (slicer[ii] == query):
                query_after = query[:-1] + chr(query[-1] + 1).encode("utf-8")
                jj = bisect_left(slicer, query_after)
                if (jj < len(self.suffix_arr) and slicer[jj] == query) or slicer[
                    jj - 1
                ] != query:
                    raise RuntimeError("index broken, should not happen")
                t2 = time.time()
                # Find row numbers
                res = (
                    numpy.searchsorted(self.offset_arr, self.suffix_arr[ii:jj] + 1) - 1
                )
                t3 = time.time()
                logger.debug(
                    "SuffixIndex.find({}) bisect: {}, zoom: {}, extract: {} --> {}rows".format(
                        query, t1 - t0, t2 - t1, t3 - t2, len(res)
                    )
                )
                return res
        return []


class SuffixIndexError(Exception):
    pass


class SuffixIndexInconsitency(SuffixIndexError):
    pass
