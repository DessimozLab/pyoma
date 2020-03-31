import tables
import numpy as np

try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache
from builtins import filter


class RelationsOfEntry(object):
    """
    This container class contains the relations from one given query
    gene in one other genome. E.g. all the BestMatches from HUMAN2 in
    MOUSE, or the StablePairs from CHICK4 in DROME... The different
    types of relations are in fact implemented as different
    subclasses.
    """

    @staticmethod
    def filter(x):
        if isinstance(x, np.ndarray):
            return np.array([True] * len(x), dtype="b")
        else:
            return np.array([True], dtype="b")

    def __init__(self, data):
        """create a RelationsOfEntry object from 'data'.
        data is assumed to be a numpy structured array containing
        matches from one query entry only, i.e. all EntryNr1 elements
        have the same value. Further, it is assumed that the data
        is ordered according to EntryNr2 (inc order)."""
        self.data = data

    def __iter__(self):
        return filter(self.filter, self.data)

    def __contains__(self, item):
        idx = np.where(self.data["EntryNr2"] == item)
        return self.filter(self.data[idx]).any()

    def relations(self):
        """return a copy of the matches restricted to the rows
        passing the filter criterion."""
        return self.data[self.filter(self.data)]

    def set_relations(self, value):
        if isinstance(value, (list, set, int)):
            value = np.array(value)
        elif not isinstance(value, np.ndarray):
            raise ValueError("Parameter not expected")

        if value.dtype == np.int:
            value_idx = self.data["EntryNr2"].searchsorted(value)
            if not (self.data["EntryNr2"][value_idx] == value).all():
                missing_matches = value[self.data["EntryNr2"][value_idx] != value]
                raise ValueError(
                    u"Not all entries found in Matches: {0!r:s}".format(missing_matches)
                )
        elif value.dtype == self.data.dtype:
            raise NotImplementedError(u"setitem with replacement view not implemented")
        else:
            raise ValueError(
                u"Unexpected dtype in parameter: {0!:s}".format(value.dtype)
            )
        self._set_relationflags_on_rows(value_idx)

    def _set_relationflags_on_rows(self, row_indexes):
        """method which sets the column flags to the appropriate
        relation type for the indexes of the data buffer. The base
        class does not change any flags. """
        pass


class CanonicalBestMatchesOfEntry(RelationsOfEntry):
    """"variant that only looks at the main splicing variant"""

    @staticmethod
    def filter(row):
        return row["isCanonicalSplicing"]

    def _set_relationflags_on_rows(self, row_indexes):
        mask = np.zeros(len(self.data), dtype="bool")
        mask[row_indexes] = True
        self.data["isCanonicalSplicing"] = mask


class StablePairsOfEntry(RelationsOfEntry):
    """Stable Pairs have to be Canonical Splicing and set isSP to true"""

    @staticmethod
    def filter(row):
        return np.logical_and(row["isSP"], row["isCanonicalSplicing"])

    def _set_relationflags_on_rows(self, row_indexes):
        mask = np.zeros(len(self.data), dtype="bool")
        mask[row_indexes] = True
        self.data["isSP"] = mask
        self.data["isCanonicalSplicing"][row_indexes] = True


class VPairsOfEntry(RelationsOfEntry):
    """Verifiy Pairs have to be Canonical Splicing and set isVP to true.
    So far we do not enforce isSP to be true as well, because of consistency
    step where we potentially augment non-sps to vps."""

    @staticmethod
    def filter(row):
        return np.logical_and(row["isVP"], row["isCanonicalSplicing"])

    def _set_relationflags_on_rows(self, row_indexes):
        mask = np.zeros(len(self.data), dtype="bool")
        mask[row_indexes] = True
        self.data["isVP"] = mask
        self.data["isCanonicalSplicing"][row_indexes] = True


class RelationManager(object):
    def __init__(self, data):
        self.data = data
        self._data_changed = False

    def __getitem__(self, key):
        if key == "SP":
            res = StablePairsOfEntry(self.data)
        elif key == "VP":
            res = VPairsOfEntry(self.data)
        elif key == "BM":
            res = CanonicalBestMatchesOfEntry(self.data)
        elif key == "ALL":
            res = RelationsOfEntry(self.data)
        else:
            raise KeyError(u"Unexpected key: {0!r:s}".format(key))
        return res

    def __setitem__(self, key, value):
        rels_obj = self.__getitem__(key)
        rels_obj.set_relations(value)
        self._data_changed = True

    def is_data_insync(self):
        return not self._data_changed


class GenomePair(object):
    def __init__(self, mtab, itab):
        self.matches = mtab
        self.entry_offset = itab
        self.rels = [None] * len(itab)

    def __getitem__(self, item):
        if not isinstance(item, int):
            raise IndexError(u"Invalid index or slice: {0!r:s}".format(item))
        if item < 0:
            raise IndexError(u"Invalid index or slice: {0!r:s}".format(item))
        if item < 0:
            raise IndexError(u"Invalid index or slice: {0!r:s}".format(item))
        if item < 0:
            item += len(self.entry_offset)
        else:
            item -= 1
        if item >= len(self.entry_offset):
            raise IndexError(u"Index is out of bound: {}".format(item))
        if not self.rels[item] is None:
            return self.rels[item]

        range_of_entry = self.entry_offset[item]
        data = self.matches[range_of_entry[0] : range_of_entry[1]]
        self.rels[item] = RelationManager(data)
        return self.rels[item]

    def flush(self):
        modified = 0
        for i, rel_man in enumerate(self.rels):
            if rel_man is not None and not rel_man.is_data_insync():
                modified += self.matches.modify_rows(
                    start=self.entry_offset[i, 0],
                    stop=self.entry_offset[i, 1],
                    rows=rel_man.data,
                )
                rel_man._data_changed = False

    def close(self):
        self.flush()
        del self.rels


class OmaDB(object):
    def __init__(self, file, mode="r"):
        if mode != "r":
            filters = tables.Filters(complevel=6, complib="zlib")
        else:
            filters = None
        self.db_handle = tables.open_file(file, mode=mode, filters=filters)
        self.genome_data = self.db_handle.get_node("/GenomeSummary")

    @lru_cache(maxsize=256)
    def matches(self, genome1, genome2):
        try:
            mtab = self.db_handle.get_node(
                "/Matches/{}/{}/Relations".format(genome1, genome2)
            )
            itab = self.db_handle.get_node(
                "/Matches/{}/{}/ProteinIndex".format(genome1, genome2)
            )
        except tables.NoSuchNodeError:
            raise KeyError(u"genome pair does not exist in database")
        gs = self.genome_data.read_where("UniProtSpeciesCode==b'{}'".format(genome1))
        if len(itab) != gs["TotEntries"]:
            raise OmaDBError(
                u"Nr of Protein does not match: {} vs {}".format(
                    len(itab), gs["TotEntries"]
                )
            )
        return GenomePair(mtab, itab)

    def close(self):
        self.db_handle.close()


class OmaDBError(Exception):
    pass
