import abc
import collections
import itertools
import logging
import re
from typing import Union, List
from dataclasses import dataclass

from .exceptions import (
    InvalidId,
    UnknownSpecies,
    InvalidTaxonId,
    AmbiguousID,
    TooUnspecificQuery,
)
from . import db, models

logger = logging.getLogger(__name__)


class BaseSearch(metaclass=abc.ABCMeta):
    PRIO = 100

    def __init__(self, pyomadb: db.Database, term: Union[int, str]):
        self.db = pyomadb
        self.term = term

    def set_max_entries(self, max_nr: int):
        pass

    def search_entries(self):
        pass

    def search_groups(self):
        pass

    def search_species(self):
        pass

    def search_ancestral_genomes(self):
        pass

    def count_entries(self):
        nr_entries = self.search_entries()
        return 0 if nr_entries is None else len(nr_entries)

    def count_groups(self):
        nr_groups = self.search_groups()
        return 0 if nr_groups is None else len(nr_groups)

    def count_species(self):
        nr_species = self.search_species()
        return 0 if nr_species is None else len(nr_species)

    def count_ancestral_genomes(self):
        nr_ancestral_genomes = self.search_ancestral_genomes()
        return 0 if nr_ancestral_genomes is None else len(nr_ancestral_genomes)


class BaseWithEntryFilter(BaseSearch):
    def __init__(self, pyomadb: db.Database, term: Union[int, str]):
        super().__init__(pyomadb, term)
        self._max_entries = None
        self._entry_nr_filter = None

    def set_max_entries(self, max_nr: int):
        self._max_entries = max_nr

    def set_entry_nr_filter(self, entry_filter: Union[set, tuple]):
        self._entry_nr_filter = entry_filter


class OmaGroupSearch(BaseSearch):
    PRIO = 50

    @models.LazyProperty
    def _matched_groups(self):
        try:
            grps = [self.db.resolve_oma_group(self.term)]
        except AmbiguousID as e:
            grps = e.candidates
        except InvalidId:
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
    PRIO = 70

    def __init__(
        self,
        pyomadb: db.Database,
        term: Union[int, str],
        level: Union[None, int] = None,
    ):
        super().__init__(pyomadb, term)
        self.outdated_query_hog = False
        self._level = None
        self._get_protein_entries = False
        self._matched_hogs = None
        if level is not None:
            try:
                self.set_taxon_filter(models.AncestralGenome(self.db, level))
            except UnknownSpecies:
                pass

    def set_taxon_filter(self, taxon):
        if not isinstance(taxon, models.AncestralGenome):
            return
        if self._level != taxon.scientific_name:
            self._level = taxon.scientific_name
            self._get_protein_entries = len(taxon.extant_genomes) < 3
            self._matched_hogs = None
            self.outdated_query_hog = False

    def _map_forward_outdated_hogid(self, hogid):
        try:
            with db.HogIdForwardMapper(self.db) as fwd:
                return fwd.map_hogid(hogid)
        except IOError as e:
            return {}

    def get_matched_hogs(self):
        if self._matched_hogs is not None:
            return self._matched_hogs

        match = re.match(
            r"(?P<id>(?P<prefix>HOG:)?(?P<rel>[A-Z]*)(?P<fam>\d+)(?P<subid>[a-z0-9.]*))(?:_(?P<taxid>\d+))?",
            self.term,
        )
        hogs = []
        if match is not None:
            level = None
            logger.debug("analyse token %s: %s", self.term, match.groupdict())
            if match.group("taxid") is not None:
                try:
                    n = self.db.tax.get_taxnode_from_name_or_taxid(match.group("taxid"))
                    level = n[0]["Name"].decode()
                except (InvalidTaxonId, KeyError):
                    pass
            new_hogs = {}
            if match.group("rel") != self.db.release_char:
                old_hogid = ("" if match.group("prefix") else "HOG:") + match.group("id")
                new_hogs = self._map_forward_outdated_hogid(old_hogid)
                self.outdated_query_hog = True
                ids = new_hogs.keys()
            else:
                ids = (
                    ("" if match.group("prefix") else "HOG:") + match.group("id"),
                    "HOG:{}{:07d}".format(self.db.release_char, int(match.group("fam"))),
                )
            if self._level is not None:
                # we have a taxon filter, search for subhogs/superhogs of ids
                ids_to_look = ids if self.outdated_query_hog else (ids[0],)
                for hog_id in ids_to_look:
                    for subhog in self.db.iter_hogs_at_level(hog_id, level=self._level):
                        hogs.append(models.HOG(self.db, subhog))
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
        self._matched_hogs = hogs
        return hogs

    def search_entries(self):
        if not self._get_protein_entries:
            return None
        return list(
            itertools.chain.from_iterable(hog.members for hog in self.get_matched_hogs() if hog.level == self._level)
        )

    def search_groups(self):
        return self.get_matched_hogs()


class GOSearch(BaseWithEntryFilter):
    PRIO = 30

    @models.LazyProperty
    def _matched_entries(self):
        try:
            if self._max_entries is None and self.estimated_occurrences > 1000:
                self.set_max_entries(1000)
            res = list(
                self.db.entrynrs_with_go_annotation(
                    self.term,
                    limit=self._max_entries,
                    entrynr_filter=self._entry_nr_filter,
                )
            )
        except InvalidId:
            res = []
        return res

    def search_entries(self):
        return [models.ProteinEntry(self.db, en) for en in self._matched_entries]

    @models.LazyProperty
    def estimated_occurrences(self):
        return self.db.count_go_annotations(self.term)

    def count_entries(self):
        return self.estimated_occurrences


class ECSearch(BaseWithEntryFilter):
    PRIO = 25

    @models.LazyProperty
    def _matched_entries(self):
        if self._max_entries is None and self.estimated_occurrences > 1000:
            self.set_max_entries(1000)
        return list(
            int(z)
            for z in self.db.entrynrs_with_ec_annotation(
                self.term, limit=self._max_entries, entrynr_filter=self._entry_nr_filter
            )
        )

    def search_entries(self):
        return [models.ProteinEntry(self.db, en) for en in self._matched_entries]

    @models.LazyProperty
    def estimated_occurrences(self):
        return self.db.count_ec_annotations(self.term)

    def count_entries(self):
        return self.estimated_occurrences


class DomainSearch(BaseWithEntryFilter):
    PRIO = 40

    @models.LazyProperty
    def _matched_entries(self):
        if self._max_entries is None and self.estimated_occurrences > 1000:
            self.set_max_entries(1000)
        return list(
            int(z)
            for z in self.db.entrynrs_with_domain_id(
                self.term, limit=self._max_entries, entrynr_filter=self._entry_nr_filter
            )
        )

    def search_entries(self):
        return [models.ProteinEntry(self.db, en) for en in self._matched_entries]

    @models.LazyProperty
    def estimated_occurrences(self):
        return self.db.count_domain_id_annotations(self.term)

    def count_entries(self):
        return self.estimated_occurrences


class TaxSearch(BaseSearch):
    PRIO = 10

    @models.LazyProperty
    def _matched_taxons(self):
        try:
            tax_nodes = self.db.tax.get_taxnode_from_name_or_taxid(self.term)
            tax = [(int(z), 1.0) for z in tax_nodes["NCBITaxonId"]]
        except KeyError:
            try:
                genome = self.db.id_mapper["OMA"].identify_genome(self.term)
                tax = [(int(genome["NCBITaxonId"]), 1.0)]
            except UnknownSpecies:
                # fuzzy match of all extant genomes in OMA
                genomes, genome_approx_scores = self.db.id_mapper["OMA"].approx_search_genomes(self.term, scores=True)
                logger.debug("%s matches approximately to %s", self.term, genomes)
                # fuzzy match of all taxonomic names in OMA.
                approx_matches = self.db.tax.approx_search(self.term)
                logger.debug("'%s' matches approximately to %s", self.term, approx_matches)
                tax_nodes = self.db.tax.get_taxnode_from_name_or_taxid([z[1] for z in approx_matches])
                tax = [(int(z["NCBITaxonId"]), approx_match[0]) for z, approx_match in zip(tax_nodes, approx_matches)]
                tax.extend((g.ncbi_taxon_id, score) for g, score in zip(genomes, genome_approx_scores))
        return tax

    def search_entries(self):
        enrs = set([])
        for g in self.search_species():
            enrs.update(range(g.entry_nr_offset + 1, g.entry_nr_offset + len(g) + 1))
        return enrs

    def search_ancestral_genomes(self):
        ancestral_genomes = []
        for tax, score in self._matched_taxons:
            ag = models.AncestralGenome(self.db, tax)
            if len(ag.extant_genomes) <= 1:
                continue
            if ag.sciname not in self.db.tax.all_hog_levels:
                orig_sciname = ag.sciname
                while True:
                    childs = self.db.tax._direct_children_taxa(int(ag.ncbi_taxon_id))
                    if len(childs) == 0 or len(childs) > 1:
                        break
                    ag = models.AncestralGenome(self.db, int(childs[0]["NCBITaxonId"]))
                ag.common_name = f"Representative for '{orig_sciname}'"
            ag.match_score = score
            ancestral_genomes.append(ag)
        return ancestral_genomes

    def search_species(self):
        extant_species = {}
        for tax, score in self._matched_taxons:
            extant_species_of_tax = models.AncestralGenome(self.db, tax).extant_genomes
            if len(extant_species_of_tax) > 1:
                score /= 2
            for sp in extant_species_of_tax:
                if sp.ncbi_taxon_id not in extant_species:
                    sp.match_score = score
                    extant_species[sp.ncbi_taxon_id] = sp
                else:
                    genome = extant_species[sp.ncbi_taxon_id]
                    genome.match_score = max(score, genome.match_score)
        return list(extant_species.values())


class SequenceSearch(BaseSearch):
    PRIO = 90

    def __init__(
        self,
        pyomadb: db.Database,
        term: str,
        strategy: Union[None, str] = None,
    ):
        super().__init__(pyomadb, term)
        self.strategy = strategy.lower() if strategy else "mixed"
        self.seq = self.db.seq_search._sanitise_seq(self.term)
        self.entry_filter = None
        self._matched_seqs = None
        self._nr_kmer_hits = None
        self._max_entries = None

    def set_entry_nr_filter(self, filter: Union[tuple, set]):
        if isinstance(filter, tuple) and len(filter) != 2:
            raise ValueError("filter parameter must be range tuple or a set")
        self.entry_filter = filter
        self._matched_seqs = None
        self._nr_kmer_hits = None

    def set_max_entries(self, max_nr: int):
        self._max_entries = max_nr

    def get_matched_seqs(self):
        if self._matched_seqs is not None:
            return self._matched_seqs

        if len(self.seq) < 5:
            logger.debug("too short sequence motif to search: %s", self.seq)
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

            self._nr_kmer_hits = len(exact_matches)

        if self.strategy == "approx" or self.strategy == "mixed" and len(res) == 0:
            approx, kmer_hits = self.db.seq_search.approx_search(
                self.seq,
                n=self._max_entries,
                is_sanitised=True,
                entrynr_range=self.entry_filter,
                return_kmer_hits=True,
            )
            for en, align_res in approx:
                pe = models.ProteinEntry(self.db, en)
                pe.mode = "approx"
                pe.score = align_res["score"]
                pe.alignment = align_res["alignment"]
                res[en] = pe

            self._nr_kmer_hits = len(kmer_hits)

        self._matched_seqs = res
        return res

    def search_entries(self):
        return list(self.get_matched_seqs().values())

    def search_groups(self):
        grps = collections.Counter(p.oma_group for p in self.get_matched_seqs().values() if p.oma_group != 0)
        return [models.OmaGroup(self.db, grp) for grp, cnt in grps.most_common(10)]

    def count_entries(self):
        if self._matched_seqs is None:
            self.get_matched_seqs()
        if self._nr_kmer_hits is None:
            return len(self._matched_seqs)
        else:
            return self._nr_kmer_hits


class XRefSearch(BaseSearch):
    PRIO = 85

    def __init__(
        self,
        pyomadb: db.Database,
        term: str,
        max_matches: Union[None, int] = None,
    ):
        super().__init__(pyomadb, term)
        self.max_matches = max_matches
        self._matched_entries = None
        self.entry_filter = None

    def set_max_entries(self, max_nr: int):
        self.max_matches = max_nr

    @models.LazyProperty
    def estimated_occurrences(self):
        return self.db.id_mapper["XRef"].search_helper.count(self.term)

    def count_entries(self):
        return self.estimated_occurrences

    def set_entry_nr_filter(self, enr_filter: Union[tuple, set]):
        if isinstance(enr_filter, tuple) and len(enr_filter) > 2:
            raise ValueError("enr_filter parameter must be range tuple or a set")
        self.entry_filter = enr_filter
        self._matched_entries = None

    def get_matched_proteins(self):
        if self._matched_entries is None:
            filt = None
            rng = None
            if self.entry_filter is not None:
                if isinstance(self.entry_filter, tuple):
                    rng = self.entry_filter
                else:
                    if len(self.entry_filter) == 0:
                        logger.info("short-circuit XRefSearch as entry_filter is empty: %s", self.entry_filter)
                        self._matched_entries = {}
                        return self._matched_entries
                    rng = (min(self.entry_filter), max(self.entry_filter))
                    filt = lambda enr: enr in self.entry_filter

            try:
                self._matched_entries = self.db.id_resolver.search_protein(
                    self.term, limit=self.max_matches, entrynr_range=rng
                )
            except TooUnspecificQuery:
                logger.exception("XRefSearch with term %s and entry_range %s", self.term, rng)
                self._matched_entries = {}
            if filt is not None:
                self._matched_entries = {enr: v for enr, v in self._matched_entries.items() if filt(enr)}
        return self._matched_entries

    def search_entries(self):
        res = []
        for enr, xrefdata in self.get_matched_proteins().items():
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
        if self.entries is not None and hasattr(other, "set_entry_nr_filter"):
            other.set_entry_nr_filter(self.entries_set)
        if (
            self.ancestral_genomes_set is not None
            and len(self.ancestral_genomes_set) == 1
            and hasattr(other, "set_taxon_filter")
        ):
            ag = next(iter(self.ancestral_genomes_set))
            other.set_taxon_filter(self.ancestral_genomes[ag])

        res = SearchResult()
        for aspect, key in zip(
            ("entries", "groups", "species", "ancestral_genomes"),
            ("entry_nr", "group_nbr", "ncbi_taxon_id", "ncbi_taxon_id"),
        ):
            term_aspect_result = getattr(other, "search_" + aspect)()
            if term_aspect_result is not None:
                if aspect == "groups" and isinstance(other, HogIDSearch):
                    key = "hog_id"
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
                    common_keys = getattr(self, aspect + "_set") & keyset
                    setattr(res, aspect + "_set", common_keys)
                    d = getattr(self, aspect)
                    if len(d) > 0:
                        common_dict = {k: d[k] for k in common_keys}
                        common_dict.update({k: keyvals[k] for k in common_keys})
                    else:
                        common_dict = keyvals
                    setattr(res, aspect, common_dict)
            else:
                setattr(res, aspect + "_set", getattr(self, aspect + "_set"))
                setattr(res, aspect, getattr(self, aspect))
        return res


def search(tokens: List[BaseSearch], entry_limit=None):
    sorted_tokens = sorted(tokens, key=lambda t: t.PRIO)
    if entry_limit is not None:
        sorted_tokens[-1].set_max_entries(entry_limit)
    res = SearchResult()
    for token in sorted_tokens:
        res &= token
    return res
