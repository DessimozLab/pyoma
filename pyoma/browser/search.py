from . import db, models
from typing import Union
import abc


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


class TaxSearch(BaseSearch):
    def fuzzy_match(self):
        pass

    @models.LazyProperty
    def _matched_taxons(self):
        try:
            tax_nodes = self.db.tax.get_taxnode_from_name_or_taxid(self.term)
            tax = [int(z) for z in tax_nodes["NCBITaxonId"]]
        except KeyError:
            try:
                extend_genome = self.db.id_mapper["OMA"].identify_genome(self.term)
                tax = [extend_genome["NCBITaxonId"]]
            except db.UnknownSpecies:
                # fuzzy match of all taxonomic names in OMA.
                tax = []
        return tax

    def search_entries(self):
        pass

    def search_species(self):
        return [models.Genome(self.db, tax) for tax in self._matched_taxons]
