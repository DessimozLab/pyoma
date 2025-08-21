import collections
import itertools
import os
import logging
from functools import lru_cache
from typing import List, Set, Optional
import re
import numpy
import tables

from ..builder import DBBuilder
from ..function_propagation_omagroup import FunctionPredictor
from ...convert import create_index_for_columns, sort_table
from ...tablefmt import GeneOntologyTable
from ...geneontology import GeneOntology, OntologyParser, AnnotationParser, GOA_Annotation, AnnotationFilter

logger = logging.getLogger(__name__)


class GeneOntologyManager:
    def __init__(self, fpath, annotation_path, ontology_path):
        self.fpath = fpath
        self.annotation_path = annotation_path
        self.ontology_path = ontology_path
        self.go = None
        self._go_buf = []
        self._inf_buf = []
        self.phase = -1  # 0 -> accept ontology, 1 -> collect annotations, 2 -> infer annotations

    def __enter__(self):
        self.h5 = tables.open_file(
            self.fpath, mode="w", filters=tables.Filters(complevel=7, complib="blosc2", fletcher32=True)
        )
        root, name = os.path.split(self.annotation_path)
        self.go_tab = self.h5.create_table(
            root, name, GeneOntologyTable, "Gene Ontology annotations", expectedrows=1e8, createparents=True
        )
        self._inf_tab = self.h5.create_table(root, name + "_infer", GeneOntologyTable, expectedrows=1e8)
        self.phase = 0
        return self

    def add_ontology(self, obo_path):
        assert self.phase == 0
        # check that ontology file is not broken. if we can build it, it should be ok
        self.go = GeneOntology(OntologyParser(obo_path))
        self.go.parse()

        with open(obo_path, "rb") as fh:
            go_obo = fh.read()
        root, name = os.path.split(self.ontology_path)
        obo = self.h5.create_carray(
            root,
            name,
            title="Gene ontology hierarchy definition",
            createparents=True,
            obj=numpy.ndarray(len(go_obo), buffer=go_obo, dtype=tables.StringAtom(1)),
        )
        obo.set_attr("ontology_release", self._get_obo_version(obo))
        self.phase = 1

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._flush_buffers()
        if self.phase == 2:
            # copy over inf_tab into go_tab
            self._bulk_copy_table(self._inf_tab, self.go_tab)
        # remove temporary _inf_tab (it has been copied over)
        self._inf_tab.remove()
        self.go_tab.flush()
        create_index_for_columns(self.go_tab, "EntryNr", "TermNr")
        sort_table(self.go_tab, col_order=["EntryNr", "TermNr", "Evidence"])

    def _get_obo_version(self, obo_arr):
        header = obo_arr[0:1000].tobytes()
        rel_info = re.search(rb"data-version:\s*(?P<version>[\w/_ -]+)", header)
        if rel_info is not None:
            rel_info = rel_info.group("version").decode()
        return rel_info

    def _flush_buffers(self):
        logger.info("flushing go annotations buffers")
        if len(self._go_buf) > 0:
            self.go_tab.append(self._go_buf)
        self._go_buf = []
        if len(self._inf_buf) > 0:
            self._inf_tab.append(self._inf_buf)
        self._inf_buf = []

    def switch_collect_to_inference_phase(self):
        """method to switch from collect phase to inference phase.

        After calling this method, you can no longer call add_annotations
        on this object. From then on, only calls to add_inference are allowed."""
        # switch to inference phase (2). build index of annotations
        assert self.phase == 1
        self._flush_buffers()
        self.go_tab.flush()
        create_index_for_columns(self.go_tab, "EntryNr", "TermNr")
        self.phase = 2

    def _bulk_copy_table(self, source_tab: tables.Table, target_tab: tables.Table, chunk: int = 1_000_000):
        # 1. remove all indexes from taget_tab
        for c in target_tab.colnames:
            if target_tab.colindexed[c]:
                target_tab.colinstances[c].remove_index()

        # 2. Copy data
        for start in range(0, source_tab.nrows, chunk):
            stop = min(start + chunk, source_tab.nrows)
            target_tab.append(source_tab[start:stop])
            logger.debug(f"copied chunk {start}:{stop} of {source_tab.nrows} rows")
        target_tab.flush()

    def annotation_generated_date(self, date):
        self.go_tab.set_attr("annotations_generated", date)

    def add_annotations(self, enrs: Set[int], anno: GOA_Annotation):
        """parse go annotations and add them to the go buffer"""
        assert self.phase == 1, "cannot add annotations if phase != 1"
        if len(enrs) == 0:
            return
        try:
            term = self.go.term_by_id(anno.term_id)
        except ValueError:
            logger.warning(f"annotation {anno}: term {anno.term_id} not valid (obsolete?)")
            return

        ev = anno.evidence.encode("utf-8")
        ref = anno.db_ref.split("|")[0].encode("utf-8")
        for enr in enrs:
            self._go_buf.append((enr, term.id, ev, ref))
        if len(self._go_buf) > 2e6:
            self._flush_buffers()

    def add_inference(self, enr: int, term_nr: int, ref="OMA_Fun:001"):
        assert self.phase == 2, "cannot add inferences if phase != 2"
        try:
            term = self.go.term_by_id(term_nr)
        except ValueError:
            logger.warning(f"annotation: term {term_nr} not valid (obsolete?)")
            return
        self._inf_buf.append((enr, term_nr, b"IEA", ref.encode("utf-8")))
        if len(self._inf_buf) > 2e6:
            self._flush_buffers()


class XRefBasedMapper:
    def __init__(self, xref_db_path: os.PathLike):
        self.xref_db_path = xref_db_path

    def __enter__(self):
        self.xref_db = tables.open_file(self.xref_db_path, mode="r")
        self.xrefs: tables.Table = self.xref_db.get_node("/XRef")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.xref_db.close()

    @lru_cache(maxsize=65536)
    def find_id(self, id_):
        if m := re.match(r"(.*)\.\d{1,3}$", id_):
            id_ = m.group(1)
        query = id_.encode("utf-8")
        condvals = {"low": query, "high": query[:-1] + bytes([query[-1] + 1])}
        stmt = "(XRefId >= low) & (XRefId < high)"
        it = self.xrefs.where(stmt, condvals)
        xrefs = numpy.fromiter(
            (row.fetch_all_fields() for row in itertools.islice(it, 100)),
            dtype=self.xrefs.dtype,
        )
        if len(xrefs) == 0:
            return set([])
        xrefs.sort(order=["Verification", "EntryNr"])
        index = numpy.searchsorted(xrefs["Verification"], xrefs[0]["Verification"], side="right")
        return set(xrefs["EntryNr"][:index])


class Stats:
    def __init__(self):
        self.counts = collections.defaultdict(lambda: [0, 0, 0, 0])

    def log(self, taxid, kind):
        if kind < 0:
            self.counts[taxid][0] += 1
        elif kind == 0:
            self.counts[taxid][1] += 1
        elif kind == 1:
            self.counts[taxid][2] += 1
        if kind > 1:
            self.counts[taxid][3] += 1

    def summary(self):
        print("Mapping statistics")
        print("=" * 80)
        print(f"{'Taxon id':>15s}|{'NOT':>13s}|{'ID not found':>15s}|{'Found 1:1 match':>17s}|{'>1 ID matched':>15s}")
        print("-" * 80)
        for taxid in sorted(self.counts.keys(), key=lambda x: -sum(self.counts[x])):
            z = self.counts[taxid]
            print(f"{taxid:15d}|{z[0]:13d}|{z[1]:15d}|{z[2]:17d}|{z[3]:15d}")
        print("-" * 80)


def import_go(
    obo: os.PathLike,
    gafs: List[os.PathLike],
    xref_db: os.PathLike,
    relevant_taxid: Set[int],
    og_db: os.PathLike,
    out: os.PathLike,
    clades: Optional[List[str]] = None,
):
    stats = Stats()
    with GeneOntologyManager(out, "/Annotations/GeneOntology", "/Ontologies/GO") as go_man:
        with XRefBasedMapper(xref_db) as mapper:
            go_man.add_ontology(obo_path=obo)
            for gaf in gafs:
                logger.info(f"importing GO annotations from {gaf}")
                parser = AnnotationParser(gaf)
                for annotation in parser:
                    # extract taxid (potentially contains also host taxon:1|taxon:1000)
                    taxid = int(annotation.taxon.split("|")[0].split(":")[1])
                    if taxid not in relevant_taxid:
                        continue
                    if AnnotationFilter.is_negated(annotation):
                        stats.log(taxid, -1)
                        continue
                    enrs = mapper.find_id(annotation.db_obj_id)
                    go_man.add_annotations(enrs, annotation)
                    stats.log(taxid, len(enrs))
                    if len(enrs) > 1:
                        logger.info(f"{annotation.db_obj_id} mapped to {len(enrs)} entries: {enrs}")

        stats.summary()
        go_man.switch_collect_to_inference_phase()

        if len(go_man.go_tab) > 0:
            logger.info("predicting GO annotations based on OMA Groups")
            with FunctionPredictor(
                db_h5_path=og_db, ontology=go_man.go, anno_tab=go_man.go_tab, clades=clades
            ) as predictor:
                for og in predictor.oma_groups:
                    for enr, go_term in predictor.annotate_group(og):
                        go_man.add_inference(enr, go_term)
        else:
            logger.info("no GO annotations found, won't predict based on OMA Groups")

    with DBBuilder(path=out, mode="append", logger=logger) as builder:
        builder.add_gene_ontology_term_cnts()
