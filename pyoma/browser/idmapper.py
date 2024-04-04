from __future__ import division, print_function, unicode_literals
import collections
import itertools
import logging
import re

import numpy
import pandas as pd
import tables
from typing import Tuple, Mapping
from .decorators import timethis
from .suffixsearch import SuffixSearcher, SuffixIndexError
from .exceptions import TooUnspecificQuery

logger = logging.getLogger(__name__)


class GeneNamesLookup:
    def __init__(self, h5: tables.File):
        tab = h5.get_node("/XRefIndex/GeneNames").read()
        self.gene_names = {x["XRefId"]: x for x in tab}
        self.lookup_tab = h5.get_node("/XRefIndex/GeneNames_lookup")

    def _low(self, item):
        if isinstance(item, str):
            item = item.encode("utf-8")
        return item.lower()

    def __contains__(self, item):
        return self._low(item) in self.gene_names

    def count(self, term):
        term = self._low(term)
        try:
            v = self.gene_names[term]
            return int(v["Length"])
        except KeyError:
            return 0

    def get_matching_xref_row_nrs(self, term, entrynr_range=None):
        term = self._low(term)
        v = self.gene_names[term]
        lookup = self.lookup_tab[v["Offset"] : v["Offset"] + v["Length"]]
        if entrynr_range is not None:
            lookup = lookup[
                numpy.where(
                    numpy.logical_and(
                        lookup["EntryNr"] >= entrynr_range[0],
                        lookup["EntryNr"] < entrynr_range[1],
                    )
                )
            ]
        return lookup["XRefRow"]


class XRefSearchHelper:
    def __init__(self, h5):
        self.gene_name_lookup = None
        self.reduced_xref_tab = None
        self._re_version = re.compile(r"(?P<base>[\w-]+)\.\d{1,2}$")
        try:
            self.gene_name_lookup = GeneNamesLookup(h5)
            self.reduced_xref_tab = h5.get_node("/XRefIndex/XRef_reduced")
        except tables.NoSuchNodeError:
            logger.warning("No reduced XRef Index and GeneName lookup found")
            pass
        self._use_reduced = self.reduced_xref_tab is not None

        self.xref_tab = h5.get_node("/XRef")
        try:
            self.fulltext_index = SuffixSearcher.from_tablecolumn(self.xref_tab, "XRefId")
        except SuffixIndexError:
            # compability mode
            idx_node = h5.get_node("/XRef_Index")
            self.fulltext_index = SuffixSearcher.from_index_node(idx_node)
        try:
            self._xref_entry_offset = h5.get_node("/XRef_Entry_offset")
        except tables.NoSuchNodeError:
            self._xref_entry_offset = None

    def _version_free_query(self, query: str):
        m = self._re_version.match(query)
        if m is not None:
            return m.group("base")
        return query

    def _query_prefix(self, query: bytes, entrynr_range=None, exact_match=False) -> Tuple[str, Mapping]:
        condvals = {}
        if exact_match:
            stmt = "(XRefId == query)"
            condvals["query"] = query
        else:
            condvals["low"] = query
            condvals["high"] = query[:-1] + bytes([query[-1] + 1])
            stmt = "(XRefId >= low) & (XRefId < high)"
        if entrynr_range is not None:
            stmt += "& (EntryNr >= enr_min) & (EntryNr < enr_max)"
            condvals.update({k: v for k, v in zip(("enr_min", "enr_max"), entrynr_range)})
        return stmt, condvals

    def _prefix_reduced_search(self, query, entrynr_range, limit=None, exact=False):
        query = self._version_free_query(query).lower().encode("utf-8")
        if query in self.gene_name_lookup:
            return self.gene_name_lookup.get_matching_xref_row_nrs(query, entrynr_range)
        red_tab_it = self.reduced_xref_tab.where(*self._query_prefix(query, entrynr_range, exact))
        xref_rows = numpy.fromiter(itertools.islice((row["XRefRow"] for row in red_tab_it), limit), dtype="i4")
        return xref_rows

    def _prefix_reducecd_count(self, query, entrynr_range=None):
        from .db import count_elements

        query = self._version_free_query(query).lower().encode("utf-8")
        cnts = self.gene_name_lookup.count(query)
        if cnts == 0:
            cnts = count_elements(self.reduced_xref_tab.where(*self._query_prefix(query, entrynr_range)))
        return cnts

    def _suffix_search(self, query, limit=None, unspecific_exception=50000):
        cnts = self.fulltext_index.count(query)
        if cnts > unspecific_exception:
            raise TooUnspecificQuery(query, cnts)
        conservative_limit = None if limit is None else 4 * limit
        return self.fulltext_index.find(query, limit=conservative_limit)

    def _suffix_count(self, query):
        return self.fulltext_index.count(query)

    def _search_reduced(self, query, mode, entrynr_range, limit):
        exact = mode == "exact"
        if mode in ("prefix", "suffix_if_no_prefix", "exact"):
            xref_rows = self._prefix_reduced_search(query, entrynr_range, limit, exact)
            xref_rows.sort()
        if mode == "suffix" or len(xref_rows) == 0 and mode == "suffix_if_no_prefix":
            xref_rows = self._search_suffix_with_constrains(query, entrynr_range, limit)
        return self.xref_tab[xref_rows]

    def _search_suffix_with_constrains(self, query, entrynr_range, limit):
        args = {}
        if limit is not None:
            unspecific_threshold = max(1000, 10 * limit) if entrynr_range is not None else 50000
            args = {"limit": limit, "unspecific_exception": unspecific_threshold}
        xref_rows = self._suffix_search(query, **args)
        xref_rows.sort()
        if entrynr_range is not None and self._xref_entry_offset is not None:
            row_low = self._xref_entry_offset[entrynr_range[0]]
            row_high = self._xref_entry_offset[entrynr_range[1] + 1]
            xref_rows = xref_rows[numpy.where(numpy.logical_and(xref_rows >= row_low, xref_rows < row_high))]
        return xref_rows

    def _search_direct(self, query, mode, entrynr_range, limit):
        query = query.encode("utf-8")
        exact = mode == "exact"
        if mode in ("prefix", "suffix_if_no_prefix", "exact"):
            print(self.xref_tab)
            print(self.xref_tab.nrows)
            it = self.xref_tab.where(*self._query_prefix(query, entrynr_range, exact_match=exact))
            xrefs = numpy.fromiter(
                (row.fetch_all_fields() for row in itertools.islice(it, limit)),
                dtype=self.xref_tab.dtype,
            )
        if mode == "suffix" or len(xrefs) == 0 and mode == "suffix_if_no_prefix":
            xref_rows = self._search_suffix_with_constrains(query, entrynr_range, limit)
            xrefs = self.xref_tab[xref_rows]
        return xrefs

    def _prefix_direct_count(self, query, entrynr_range=None):
        from .db import count_elements

        query = query.encode("utf-8")
        it = self.xref_tab.where(*self._query_prefix(query, entrynr_range))
        return count_elements(it)

    def search(self, query, mode=None, entrynr_range=None, limit=None):
        if mode is None:
            mode = "suffix_if_no_prefix"
        if mode not in ("suffix", "prefix", "suffix_if_no_prefix", "exact"):
            raise ValueError("invalid search mode: {}".format(mode))
        if self._use_reduced:
            return self._search_reduced(query, mode, entrynr_range, limit)
        else:
            return self._search_direct(query, mode, entrynr_range, limit)

    def count(self, query, mode=None):
        cnts = 0
        if self._use_reduced:
            if mode is None or mode in ("prefix", "suffix_if_no_prefix"):
                cnts = self._prefix_reducecd_count(query)
            if mode == "suffix" or cnts == 0 and mode in (None, "suffix_if_no_prefix"):
                cnts = self._suffix_count(query)
        else:
            if mode is None or mode in ("prefix", "suffix_if_no_prefix"):
                cnts = self._prefix_direct_count(query)
            if mode == "suffix" or cnts == 0 and mode in (None, "suffix_if_no_prefix"):
                cnts = self._suffix_count(query)
        return cnts


class NoSearchXrefIdMapper(object):
    def __init__(self, db, sources=None):
        self._db = db
        self.xref_tab = db.get_hdf5_handle().get_node("/XRef")
        self.xrefEnum = self.xref_tab.get_enum("XRefSource")
        self.idtype = frozenset(nr for src, nr in self.xrefEnum._names.items() if sources is None or src in sources)
        self.verif_enum = self.xref_tab.get_enum("Verification")
        self._max_verif_for_mapping_entrynrs = 1000  # allow all verification values
        try:
            self._xref_entry_offset = db.get_hdf5_handle().get_node("/XRef_EntryNr_offset")
        except tables.NoSuchNodeError:
            self._xref_entry_offset = None

    def map_entry_nr(self, entry_nr):
        """returns the XRef entries associated with the query protein.

        The types of XRefs that are returned depends on the idtype
        class member variable. In the base-class, idtype contains
        all valid xref types. Typically, subclasses of XrefIdMapper
        will change this set.

        :param entry_nr: the numeric id of the query protein.
        :returns: list of dicts with 'source' and 'xref' keys."""
        it = self._iter_xref_for_entry_nr(entry_nr)
        res = [
            {
                "source": self.xrefEnum(row["XRefSource"]),
                "xref": row["XRefId"].decode(),
                "seq_match": self.verif_enum(row["Verification"]),
            }
            for row in it
            if row["XRefSource"] in self.idtype
        ]
        return res

    def _iter_xref_for_entry_nr(self, entry_nr):
        if self._xref_entry_offset is not None:
            start, stop = (self._xref_entry_offset[z] for z in (entry_nr, entry_nr + 1))
            it = self.xref_tab.where(
                "Verification <= {:d}".format(self._max_verif_for_mapping_entrynrs),
                start=start,
                stop=stop,
            )
        else:
            it = self.xref_tab.where(
                "(EntryNr=={:d}) & (Verification <= {:d})".format(entry_nr, self._max_verif_for_mapping_entrynrs)
            )
        return it

    def canonical_source_order(self):
        """returns the list of xref sources in order of their importance.

        Most important source - in the base class for example UniProtKB/SwissProt
        are first. The canonical order is defined in the enum definition.

        :returns: list of source strings"""
        return [self.xrefEnum(z) for z in sorted(self.idtype)]

    def iter_xrefs_for_entry_nr(self, entry_nr):
        """Iterate over the xrefs of a given entry number.

        This method returns a dict with 'source' and 'xref' fields
        (both str) holding the information of the xref record.

        :param entry_nr: the numeric id of the query protein"""
        for row in self._iter_xref_for_entry_nr(entry_nr):
            if row["XRefSource"] in self.idtype:
                yield {
                    "source": self.xrefEnum(row["XRefSource"]),
                    "xref": row["XRefId"].decode(),
                }

    def _combine_query_values(self, field, values):
        parts = ["({}=={})".format(field, z) for z in values]
        return "(" + "|".join(parts) + ")"

    def map_many_entry_nrs(self, entry_nrs):
        """map several entry_nrs with as few db queries as possible
        to their cross-references. The function returns a
        :class:`numpy.recarray` containing all fields as defined in
        the table.

        :param entry_nrs: a list with numeric protein entry ids"""
        mapped_junks = []
        chunk_size = 32
        source_condition = None
        verif_condition = None
        if len(self.idtype) < len(self.xrefEnum):
            chunk_size -= len(self.idtype)  # respect max number of condition variables.
            source_condition = self._combine_query_values("XRefSource", self.idtype)
        if self._max_verif_for_mapping_entrynrs < max(x[1] for x in self.verif_enum):
            chunk_size -= 1
            verif_condition = "(Verification <= {:d})".format(self._max_verif_for_mapping_entrynrs)
        for start in range(0, len(entry_nrs), chunk_size):
            condition_list = [self._combine_query_values("EntryNr", entry_nrs[start : start + chunk_size])]
            if source_condition:
                condition_list.append(source_condition)
            if verif_condition:
                condition_list.append(verif_condition)
            condition = " & ".join(condition_list)
            mapped_junks.append(self.xref_tab.read_where(condition))
        return numpy.lib.recfunctions.stack_arrays(mapped_junks, usemask=False)

    def map_entry_nr_range(self, start, stop):
        """maps for a whole range the entry numbers to the xrefs

        param start: the first entry nr to be mapped
        type start: int, numpy.int
        param stop: the first entry nr that is not included in the result (exlusive)
        type stop: int, numpy.int

        returns: all the mapped xrefs that are respecting
                 the filtering conditions of the actual subtype of
                 XRefIdMapper.
        rtype: :class:`numpy.lib.recarray`"""
        conditions = ["(EntryNr >= {:d}) & (EntryNr < {:d})".format(start, stop)]
        if len(self.idtype) < len(self.xrefEnum):
            source_condition = self._combine_query_values("XRefSource", self.idtype)
            conditions.append(source_condition)
        if self._max_verif_for_mapping_entrynrs < max(x[1] for x in self.verif_enum):
            verif_condition = "(Verification <= {:d})".format(self._max_verif_for_mapping_entrynrs)
            conditions.append(verif_condition)
        query = " & ".join(conditions)
        it = self.xref_tab.where(query)
        try:
            first = next(it)
        except StopIteration:
            return numpy.array([], dtype=self.xref_tab.dtype)
        res = numpy.fromiter(
            map(lambda row: row.fetch_all_fields(), itertools.chain([first], it)),
            dtype=self.xref_tab.dtype,
        )
        return res

    def source_as_string(self, source):
        """string representation of xref source enum value

        this auxiliary method converts the numeric value of
        a xref source into a string representation.

        :param int source: numeric value of xref source"""
        return self.xrefEnum(source)

    def verification_as_string(self, verif):
        """string representation of xref verifiction enum value"""
        return self.verif_enum(verif)

    def xreftab_to_dict(self, tab):
        """convert a xreftable to a dictionary per entry_nr.

        All rows in `tab` are converted into a nested dictionary
        where the outer key is a protein entry number and the
        inner key the xref source type.

        :param tab: a :class:`numpy.recarray` corresponding to XRef
            table definition to be converted"""
        xrefdict = collections.defaultdict(dict)
        for row in tab:
            try:
                typ = self.xrefEnum(row["XRefSource"])
            except IndexError:
                logger.warning("invalid XRefSource value in %s", row)
                continue
            if typ not in xrefdict[row["EntryNr"]]:
                xrefdict[row["EntryNr"]][typ] = {
                    "id": row["XRefId"],
                    "seq_match": self.verif_enum(row["Verification"]),
                }
        return xrefdict


class XrefIdMapper(NoSearchXrefIdMapper):
    def __init__(self, db):
        super(XrefIdMapper, self).__init__(db)
        self.search_helper = XRefSearchHelper(db.get_hdf5_handle())

    @timethis(logging.DEBUG)
    def search_xref(self, xref, is_prefix=False, match_any_substring=False):
        """identify proteins associcated with `xref`.

        The crossreferences are limited to the types in the class
        member `idtype`. In the base class, all types are valid
        xrefs. The method returns a :class:`numpy.recarry` defined
        for the XRef table with all entries pointing to `xref`.

        The method by default returns only exact matches. By setting
        `is_prefix` to True, one can indicated that the requested xref
        should be interpreted as a prefix and all entries matching this
        prefix should be returned.

        :param str xref: an xref to be located
        :param bool is_prefix: treat xref as a prefix and return
                     potentially several matching xrefs
        :param bool match_any_substring: use a suffix index to find
                     any occurrence of the query xref (not limited
                     to the beginning)"""

        if match_any_substring:
            mode = "suffix"
        elif is_prefix:
            mode = "prefix"
        else:
            mode = "exact"
        res = self.search_helper.search(xref, mode=mode)
        if len(res) > 0 and len(self.idtype) < len(self.xrefEnum):
            res = res[numpy.in1d(res["XRefSource"], list(self.idtype))]
        return res

    @timethis(logging.INFO)
    def search_id(self, query, limit=None, entrynr_range=None):
        source_filter = None
        try:
            prefix, term = query.split(":", maxsplit=1)
            if prefix in self.xrefEnum:
                source_filter = self.xrefEnum[prefix]
            else:
                term = query
        except ValueError:
            term = query

        result = collections.defaultdict(dict)
        for xref in self.search_helper.search(term, entrynr_range=entrynr_range, limit=limit):
            if entrynr_range is not None and not (entrynr_range[0] <= xref["EntryNr"] <= entrynr_range[1]):
                continue
            if not source_filter or xref["XRefSource"] == source_filter:
                source = self.xrefEnum(xref["XRefSource"])
                try:
                    result[xref["EntryNr"]][source].append(xref["XRefId"].decode())
                except KeyError:
                    result[xref["EntryNr"]][source] = [xref["XRefId"].decode()]
            if limit is not None and len(result) >= limit:
                break
        return result


class XRefNoApproximateIdMapper(XrefIdMapper):
    def __init__(self, db):
        super(XRefNoApproximateIdMapper, self).__init__(db)
        self._max_verif_for_mapping_entrynrs = self.verif_enum["unchecked"]


class UniProtIdMapper(XrefIdMapper):
    def __init__(self, db):
        super(UniProtIdMapper, self).__init__(db)
        self.idtype = frozenset([self.xrefEnum[z] for z in ["UniProtKB/SwissProt", "UniProtKB/TrEMBL"]])


class LinkoutIdMapper(XrefIdMapper):
    def __init__(self, db):
        super(LinkoutIdMapper, self).__init__(db)
        self.idtype = frozenset(
            [
                self.xrefEnum[z]
                for z in [
                    "UniProtKB/SwissProt",
                    "UniProtKB/TrEMBL",
                    "Ensembl Protein",
                    "Ensembl Gene",
                    "EntrezGene",
                ]
            ]
        )

    def url(self, typ, id_):
        # TODO: improve url generator in external module with all xrefs
        url = None
        try:
            id_ = id_.decode()
        except AttributeError:
            pass

        if typ.startswith("UniProtKB"):
            url = "http://uniprot.org/uniprot/{}".format(id_)
        elif typ == "EntrezGene":
            url = "http://www.ncbi.nlm.nih.gov/gene/{}".format(id_)
        elif typ.startswith("Ensembl"):
            url = "http://ensembl.org/id/{}".format(id_)
        return url

    def xreftab_to_dict(self, tab):
        xref = super(LinkoutIdMapper, self).xreftab_to_dict(tab)
        for d in list(xref.values()):
            for typ, elem in list(d.items()):
                elem["url"] = self.url(typ, elem["id"])
        return xref

    def iter_xrefs_for_entry_nr(self, entry_nr):
        """same as base clase but includes also the url as a field"""
        for xref in super(LinkoutIdMapper, self).iter_xrefs_for_entry_nr(entry_nr):
            xref["url"] = self.url(xref["source"], xref["xref"])
            yield xref


class GeneNameOrSymbolIdMapper(XRefNoApproximateIdMapper):
    def __init__(self, db):
        super(GeneNameOrSymbolIdMapper, self).__init__(db)
        self.order = [
            "Gene Name",
            "UniProtKB/SwissProt",
            "UniProtKB/TrEMBL",
            "HGNC",
            "SourceID",
        ]
        self.idtype = frozenset(self.xrefEnum[z] for z in self.order)

    def canonical_source_order(self):
        return self.order


class GeneNameAndMainDbIdMapper(XRefNoApproximateIdMapper):
    def __init__(self, db):
        super(GeneNameAndMainDbIdMapper, self).__init__(db)
        self.order = [
            "Gene Name",
            "UniProtKB/SwissProt",
            "UniProtKB/TrEMBL",
            "Ensembl Gene",
            "RefSeq",
            "EntrezGene",
            "SourceID",
        ]
        self.idtype = frozenset(self.xrefEnum[z] for z in self.order)

    def canonical_source_order(self):
        return self.order


class DomainNameIdMapper(object):
    def __init__(self, db):
        self.domain_src = db.get_hdf5_handle().root.Annotations.DomainDescription.read()
        self.domain_src.sort(order="DomainId")

    def _get_dominfo(self, domain_id):
        idx = self.domain_src["DomainId"].searchsorted(domain_id)
        if self.domain_src[idx]["DomainId"] != domain_id:
            raise KeyError("no domain info available for {}".format(domain_id))
        return self.domain_src[idx]

    def get_info_dict_from_domainid(self, domain_id):
        info = self._get_dominfo(domain_id)
        return {
            "name": info["Description"].decode(),
            "source": info["Source"].decode(),
            "domainid": domain_id.decode(),
        }
