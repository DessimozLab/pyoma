import collections
import itertools
import os
import logging
from functools import lru_cache
from typing import List, Set
import re
import numpy
import tables
from ...tablefmt import GeneOntologyTable
from ...geneontology import GeneOntology, OntologyParser, AnnotationParser, GOA_Annotation, AnnotationFilter

logger = logging.getLogger(__name__)


class GeneOntologyManager:
    def __init__(self, fpath, annotation_path, ontology_path):
        self.fpath = fpath
        self.annotation_path = annotation_path
        self.ontology_path = ontology_path
        self._go_buf = []

    def __enter__(self):
        self.h5 = tables.open_file(
            self.fpath, mode="w", filters=tables.Filters(complevel=5, complib="blosc2", fletcher32=True)
        )
        root, name = os.path.split(self.annotation_path)
        self.go_tab = self.h5.create_table(
            root, name, GeneOntologyTable, "Gene Ontology annotations", expectedrows=1e8, createparents=True
        )
        return self

    def add_ontology(self, obo_path):
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

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._flush_buffers()
        self.go_tab.flush()

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

    def annotation_generated_date(self, date):
        self.go_tab.set_attr("annotations_generated", date)

    def add_annotations(self, enrs: Set[int], anno: GOA_Annotation):
        """parse go annotations and add them to the go buffer"""
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


class XRefBasedMapper:
    def __init__(self, xref_db_path: os.PathLike):
        self.xref_db_path = xref_db_path

    def __enter__(self):
        self.xref_db = tables.open_file(self.xref_db_path, mode="r")
        self.xrefs = self.xref_db.get_node("/XRef")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.xref_db.close()

    @lru_cache(maxsize=65536)
    def find_id(self, id_):
        if m := re.match(r"(.*)\.\d{1,3}$", id):
            id_ = m.group(1)
        query = id_.encode("utf-8")
        condvals = {"low": query, "high": query[:-1] + bytes([query[-1] + 1])}
        stmt = "(XRefId >= low) & (XRefId < high)"
        it = self.xrefs.where(stmt, condvals)
        xrefs = numpy.fromiter(
            (row.fetch_all_fields() for row in itertools.islice(it, 100)),
            dtype=self.xrefs.dtype,
        )
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
            logger.info("mapped to more than one enr")

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
    obo: os.PathLike, gafs: List[os.PathLike], xref_db: os.PathLike, relevant_taxid: Set[int], out: os.PathLike
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
                    for enr in enrs:
                        go_man.add_annotations(enr, annotation)
                    stats.log(taxid, len(enrs))
    stats.summary()
