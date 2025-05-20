from logging import getLogger
from typing import Optional, List, Dict, Tuple, Set, Iterable, Union
import os
import collections

import numpy
from numpy.typing import NDArray
import tables

from ..db import Taxonomy
from ..geneontology import GeneOntology, GOterm, AnnotationFilter


logger = getLogger(__name__)
CLADES = [
    "Bacteria",
    "Archaea",
    "Fungi",
    "Viridiplantae",
    "Nematoda",
    "Arthropoda",
    "Mammalia",
    "Sauria",
    "Dictyostelium",
    "Clupeocephala",
    "Amphibia",
]


def map_clades_to_entrynr_ranges(h5: tables.File, clades: List[str]) -> NDArray:
    clade2range = {}
    tax = Taxonomy(h5.get_node("/Taxonomy").read())
    genomes = h5.get_node("/Genome").read()
    for clade in clades:
        try:
            node = tax.get_taxnode_from_name_or_taxid(clade)
        except KeyError:
            logger.warning(f"Clade {clade} not found. Will not collect any GO terms for this clade.")
            continue

        subtax = tax.get_subtaxonomy_rooted_at(node[0]["NCBITaxonId"])
        taxids = subtax.get_taxid_of_extent_genomes()
        genomes_of_clade = genomes[numpy.isin(genomes["NCBITaxonId"], taxids)]
        k = genomes_of_clade["EntryOff"].argmax()
        rng = (genomes_of_clade["EntryOff"].min(), genomes_of_clade[k]["EntryOff"] + genomes_of_clade[k]["TotEntries"])
        clade2range[clade] = rng
    rng_tab = numpy.array(
        [(rng[0], rng[1], clade) for clade, rng in clade2range.items()],
        dtype=[("Low", "i4"), ("High", "i4"), ("Clade", "U255")],
    )
    rng_tab.sort(order="Low")
    return rng_tab


def collect_terms_per_clade(go_anno_tab: tables.Table, ontology: GeneOntology, clade_ranges: NDArray):
    clade2terms = {}
    trusted_evidence = frozenset(
        [x.encode("utf-8") for x in AnnotationFilter.EXP_CODES.union(AnnotationFilter.PHYL_CODES)]
    )
    for clade in clade_ranges:
        terms = set([])
        query = "(EntryNr > low) & (EntryNr <= high)"
        cond = {"low": clade["Low"], "high": clade["High"]}
        for an in go_anno_tab.where(query, cond):
            if an["Evidence"] in trusted_evidence:
                terms.add(an["TermNr"])

        # add implicit parent terms of all annotations
        superterms = set([])
        for term in terms:
            superterms.update(ontology.get_superterms_incl_queryterm(term))
        clade2terms[clade["Clade"]] = superterms
        logger.info(f"collected {len(terms)} ({len(superterms)} including superterms) in {clade['Clade']}")
    return clade2terms


class FunctionPredictor:
    def __init__(
        self,
        db_h5_path: os.PathLike,
        anno_tab: tables.Table,
        ontology: GeneOntology,
        clades: Optional[List[str]] = None,
    ):
        self.db_h5_path = db_h5_path
        self.anno_tab = anno_tab
        self.ontology = ontology
        self.clades = clades if clades is not None else CLADES

    def __enter__(self):
        with tables.open_file(self.db_h5_path) as db:
            self.clade_ranges = map_clades_to_entrynr_ranges(db, self.clades)
            self.oma_groups = self.load_oma_groups(db.get_node("/Protein/Entries"))
        self.clade2terms = collect_terms_per_clade(self.anno_tab, self.ontology, self.clade_ranges)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def load_oma_groups(self, tab: tables.Table) -> Dict[int, List[int]]:
        grps = collections.defaultdict(list)
        for row in tab.iterrows():
            if row["OmaGroup"] > 0:
                grps[row["OmaGroup"]].append(row["EntryNr"])
        return grps

    def enr2clade(self, enrs: Union[int, NDArray[int]]) -> NDArray[numpy.str_]:
        if self.clade_ranges.size > 0:
            k = self.clade_ranges["Low"].searchsorted(enrs)
            return numpy.where(
                (k > 0) & (self.clade_ranges[k - 1]["High"] >= enrs), self.clade_ranges[k - 1]["Clade"], ""
            )
        else:
            return numpy.zeros_like(enrs, numpy.str_)

    def reliable_annotations_of_entry(self, enr: int) -> Set[GOterm]:
        terms = set()
        for row in self.anno_tab.where("EntryNr == enr", {"enr": enr}):
            if row["Evidence"] != b"IEA" or row["Reference"].decode() in AnnotationFilter.TRUST_IEA_REFS:
                terms.add(self.ontology.term_by_id(row["TermNr"]))
        return terms

    def implied_terms(self, terms: Iterable[Union[int, GOterm]]) -> Set[GOterm]:
        # add implicit parent terms of all annotations
        super_terms = set([])
        for term in terms:
            super_terms.update(self.ontology.get_superterms_incl_queryterm(term))
        return super_terms

    def annotate_group(self, grp):
        # load annotations of group members
        annos = collections.defaultdict(int)
        entry_anno = {}
        for enr in self.oma_groups[grp]:
            entry_anno[enr] = self.reliable_annotations_of_entry(enr)
            for term in self.implied_terms(entry_anno[enr]):
                annos[term] += 1

        # filter common set for sufficient prevalence
        annos = {term for term, cnt in annos.items() if cnt > min(3, len(self.oma_groups[grp]))}
        logger.debug(f"found {len(annos)} relevant GO Terms: {[str(z) for z in annos]}")

        # take the most specific annotations only
        specific = {a for a in annos if self.ontology.get_subterms(a, include_query=False).isdisjoint(annos)}
        logger.debug(f"found {len(annos)} specific GO Terms: {[str(z) for z in annos]}")

        # annotate all group members
        clades = self.enr2clade(list(entry_anno.keys()))
        logger.debug(f"clades: {clades}")
        for (enr, annos_entry), clade in zip(entry_anno.items(), clades):
            if not clade:
                continue
            candidates = specific - annos_entry
            logger.debug(f"enr {enr} [{clade}] candidates: {[str(t) for t in candidates]}")
            for candidate in candidates:
                if candidate in self.clade2terms[clade]:
                    yield enr, candidate.id
                else:
                    logger.debug(f"{candidate} not in clade terms.")
                    # if the candidate itself is not among the clade terms, check its super terms.
                    possible_implied = self.ontology.get_superterms_incl_queryterm(candidate).intersection(
                        self.clade2terms[clade]
                    )
                    specific_implied = {
                        a
                        for a in possible_implied
                        if self.ontology.get_subterms(a, include_query=False).isdisjoint(possible_implied)
                    }
                    for term in specific_implied - annos_entry:
                        yield enr, term.id
