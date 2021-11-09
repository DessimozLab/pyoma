from . import db, models
from typing import Union
import abc
import collections
import logging

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
    @models.LazyProperty
    def _matched_hogs(self):
        pass


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

    def search_ancestral(self):
        return [models.AncestralGenome(self.db, tax) for tax in self._matched_taxons]

    def search_species(self):
        sp_set = set([])
        for tax in self._matched_taxons:
            sp_set |= set(models.AncestralGenome(self.db, tax).extant_genomes)
        return sp_set


class SequenceSearch(BaseSearch):
    def __init__(
        self,
        pyomadb: db.Database,
        term: Union[int, str],
        strategy: Union[None, str] = None,
    ):
        super().__init__(pyomadb, term)
        self.strategy = strategy.lower() if strategy else "mixed"
        self.seq = self.db.seq_search._sanitise_seq(self.term)
        self.entry_filter = None

    def set_entry_nr_filter(self, filter: Union[tuple, set]):
        if isinstance(filter, tuple) and len(filter) > 2:
            raise ValueError("filter parameter must be range tuple or a set")
        self.entry_filter = filter

    @models.LazyProperty
    def _matched_seqs(self):
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
        return res

    def search_entries(self):
        return list(self._matched_seqs.values())

    def search_groups(self):
        grps = collections.Counter(
            p.oma_group for p in self._matched_seqs.values() if p.oma_group != 0
        )
        return [models.OmaGroup(self.db, grp) for grp, cnt in grps.most_common(10)]


class XRefSearch(BaseSearch):
    def __init__(
        self,
        pyomadb: db.Database,
        term: Union[int, str],
        max_matches: Union[None, int] = None,
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
