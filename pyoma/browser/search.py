import abc
import collections
import itertools
import logging
import re
from typing import Union
from dataclasses import dataclass
from . import db, models

logger = logging.getLogger(__name__)


class BaseSearch(metaclass=abc.ABCMeta):
    def __init__(self, pyomadb: db.Database, term: Union[int, str]):
        self.db = pyomadb
        self.term = term

    def search_entries(self):
        pass

    def search_groups(self):
        pass

    def search_species(self):
        pass

    def search_ancestral_genomes(self):
        pass


class OmaGroupSearch(BaseSearch):
    @models.LazyProperty
    def _matched_groups(self):
        try:
            grps = [self.db.resolve_oma_group(self.term)]
        except db.AmbiguousID as e:
            grps = e.candidates
        except db.InvalidId:
            grps = []
        return grps

    def search_groups(self):
        return [models.OmaGroup(self.db, gr) for gr in self._matched_groups]

    def search_entries(self):
        res = []
        for grp in self._matched_groups:
            for en in self.db.oma_group_members(grp):
                res.append(models.ProteinEntry(self.db, en))
        return res


class HogIDSearch(BaseSearch):
    def __init__(self, pyomadb: db.Database, term: Union[int, str]):
        super().__init__(pyomadb, term)
        self.outdated_query_hog = False

    def _map_forward_outdated_hogid(self, hogid):
        try:
            with db.HogIdForwardMapper(self.db) as fwd:
                return fwd.map_hogid(hogid)
        except IOError as e:
            return {}

    @models.LazyProperty
    def _matched_hogs(self):
        match = re.match(
            r"(?P<id>(?P<prefix>HOG:)?(?P<rel>[A-Z]*)(?P<fam>\d+)(?P<subid>[a-z0-9.]*))(?:_(?P<taxid>\d+))?",
            self.term,
        )
        hogs = []
        if match is not None:
            level = None
            if match.group("taxid") is not None:
                try:
                    n = self.db.tax.get_taxnode_from_name_or_taxid(match.group("taxid"))
                    level = n[0]["Name"].decode()
                except (db.InvalidTaxonId, KeyError):
                    pass
            new_hogs = {}
            if match.group("rel") != self.db.release_char and match.group("prefix"):
                new_hogs = self._map_forward_outdated_hogid(match.group("id"))
                self.outdated_query_hog = True
                ids = new_hogs.keys()
            else:
                ids = (
                    match.group("id"),
                    "HOG:{}{:07d}".format(
                        self.db.release_char, int(match.group("fam"))
                    ),
                )
            levels = [level]
            if level is not None:
                levels.append(None)
            for hog_id, lev in itertools.product(ids, levels):
                try:
                    hog = models.HOG(self.db, self.db.get_hog(hog_id, lev))
                    hogs.append(hog)
                    if self.outdated_query_hog:
                        hog.query_jaccard_similarity = new_hogs[hog_id]
                    else:
                        break  # not an outdated hog. break at the best candidate
                except ValueError:
                    pass
        return hogs

    def search_entries(self):
        if (
            len(self._matched_hogs) == 0
            or max(h.nr_member_genes for h in self._matched_hogs) > 100
        ):
            return None
        return list(
            itertools.chain.from_iterable(hog.members for hog in self._matched_hogs)
        )

    def search_groups(self):
        return self._matched_hogs


class GOSearch(BaseSearch):
    @models.LazyProperty
    def _matched_entries(self):
        try:
            res = list(self.db.entrynrs_with_go_annotation(self.term))
        except db.InvalidId:
            res = []
        return res

    def search_entries(self):
        return [models.ProteinEntry(self.db, en) for en in self._matched_entries]


class ECSearch(BaseSearch):
    @models.LazyProperty
    def _matched_entries(self):
        return list(int(z) for z in self.db.entrynrs_with_ec_annotation(self.term))

    def search_entries(self):
        return [models.ProteinEntry(self.db, en) for en in self._matched_entries]


class TaxSearch(BaseSearch):
    @models.LazyProperty
    def _matched_taxons(self):
        try:
            tax_nodes = self.db.tax.get_taxnode_from_name_or_taxid(self.term)
            tax = [int(z) for z in tax_nodes["NCBITaxonId"]]
        except KeyError:
            try:
                genome = self.db.id_mapper["OMA"].identify_genome(self.term)
                tax = [int(genome["NCBITaxonId"])]
            except db.UnknownSpecies:
                # fuzzy match of all taxonomic names in OMA.
                approx_matches = self.db.tax.approx_search(self.term)
                logger.debug(
                    "'{}' matches approximately to {}".format(self.term, approx_matches)
                )
                tax_nodes = self.db.tax.get_taxnode_from_name_or_taxid(
                    [z[1] for z in approx_matches]
                )
                tax = [int(z["NCBITaxonId"]) for z in tax_nodes]
        return tax

    def search_entries(self):
        enrs = set([])
        for g in self.search_species():
            enrs.update(range(g.entry_nr_offset + 1, g.entry_nr_offset + len(g) + 1))
        return enrs

    def search_ancestral_genomes(self):
        return [models.AncestralGenome(self.db, tax) for tax in self._matched_taxons]

    def search_species(self):
        sp_set = set([])
        for tax in self._matched_taxons:
            sp_set |= set(models.AncestralGenome(self.db, tax).extant_genomes)
        return list(sp_set)


class SequenceSearch(BaseSearch):
    def __init__(
        self, pyomadb: db.Database, term: str, strategy: Union[None, str] = None,
    ):
        super().__init__(pyomadb, term)
        self.strategy = strategy.lower() if strategy else "mixed"
        self.seq = self.db.seq_search._sanitise_seq(self.term)
        self.entry_filter = None
        self._matched_seqs = None

    def set_entry_nr_filter(self, filter: Union[tuple, set]):
        if isinstance(filter, tuple) and len(filter) > 2:
            raise ValueError("filter parameter must be range tuple or a set")
        self.entry_filter = filter
        self._matched_seqs = None

    def get_matched_seqs(self):
        if self._matched_seqs is not None:
            return self._matched_seqs

        if len(self.seq) < 5:
            logger.debug("too short sequence motif to search: {}".format(self.seq))
            raise ValueError("too short sequence motif")
        if self.strategy not in ("exact", "approx", "mixed"):
            raise ValueError("Invalid search strategy parameter")

        res = {}
        if self.strategy in ("exact", "mixed"):
            exact_matches = self.db.seq_search.exact_search(
                self.seq,
                only_full_length=False,
                is_sanitised=True,
                entrynr_range=self.entry_filter,
            )
            for en in exact_matches:
                pe = models.ProteinEntry(self.db, en)
                pe.mode = "exact"
                pe.score = 5 * len(self.seq)
                pe.alignment = self.seq
                res[en] = pe
        if self.strategy == "approx" or self.strategy == "mixed" and len(res) == 0:
            approx = self.db.seq_search.approx_search(
                self.seq, is_sanitised=True, entrynr_range=self.entry_filter
            )
            for en, align_res in approx:
                pe = models.ProteinEntry(self.db, en)
                pe.mode = "approx"
                pe.score = align_res["score"]
                pe.alignment = align_res["alignment"]
                res[en] = pe
        self._matched_seqs = res
        return res

    def search_entries(self):
        return list(self.get_matched_seqs().values())

    def search_groups(self):
        grps = collections.Counter(
            p.oma_group for p in self.get_matched_seqs().values() if p.oma_group != 0
        )
        return [models.OmaGroup(self.db, grp) for grp, cnt in grps.most_common(10)]


class XRefSearch(BaseSearch):
    def __init__(
        self, pyomadb: db.Database, term: str, max_matches: Union[None, int] = None,
    ):
        super().__init__(pyomadb, term)
        self.max_matches = max_matches

    @models.LazyProperty
    def estimated_occurrences(self):
        return self.db.id_mapper["XRef"].xref_index.count(self.term)

    @models.LazyProperty
    def _matched_proteins(self):
        return self.db.id_resolver.search_protein(self.term, limit=self.max_matches)

    def search_entries(self):
        res = []
        for enr, xrefdata in self._matched_proteins.items():
            p = models.ProteinEntry(self.db, enr)
            p.xref_data = xrefdata
            res.append(p)
        return res


@dataclass
class SearchResult:
    entries: dict = None
    groups: dict = None
    species: dict = None
    ancestral_genomes: dict = None
    entries_set: set = None
    groups_set: set = None
    species_set: set = None
    ancestral_genomes_set: set = None

    def __and__(self, other: BaseSearch):
        assert isinstance(other, BaseSearch)
        if hasattr(other, "set_entry_nr_filter") and self.entries is not None:
            other.set_entry_nr_filter(self.entries_set)

        res = SearchResult()
        for aspect, key in zip(
            ("entries", "groups", "species", "ancestral_genomes"),
            ("entry_nr", "group_nbr", "ncbi_taxon_id", "ncbi_taxon_id"),
        ):
            term_aspect_result = getattr(other, "search_" + aspect)()
            if term_aspect_result is not None:
                if isinstance(term_aspect_result, set):
                    keyset = term_aspect_result
                    keyvals = collections.defaultdict(dict)
                else:
                    keyset = set(getattr(z, key) for z in term_aspect_result)
                    keyvals = {getattr(z, key): z for z in term_aspect_result}

                if getattr(self, aspect) is None:
                    setattr(res, aspect, keyvals)
                    setattr(res, aspect + "_set", keyset)
                else:
                    setattr(
                        res, aspect + "_set", getattr(self, aspect + "_set") & keyset
                    )
                    setattr(res, aspect, getattr(self, aspect) | keyvals)
            else:
                setattr(res, aspect + "_set", getattr(self, aspect + "_set"))
                setattr(res, aspect, getattr(self, aspect))
        return res
