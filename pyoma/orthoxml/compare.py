from xml.etree.ElementTree import XMLParser
from datasketch import MinHash, MinHashLSH
from ..common import auto_open
import logging

logger = logging.getLogger(__name__)


class Comparer(object):
    def __init__(self):
        pass


class GeneMapper(object):
    attr = "protId"

    def map(self, attributes):
        return attributes[self.attr]


class FilteredSpeciesGeneMapper(GeneMapper):
    def __init__(self, skip_prefixes):
        prefix_lens = set(len(z) for z in skip_prefixes)
        if len(prefix_lens) != 1:
            raise ValueError("prefixes should be all equal size")
        self.skip_prefixes = skip_prefixes
        self.prefix_len = len(skip_prefixes[0])

    def map(self, attributes):
        if attributes[self.attr][: self.prefix_len] in self.skip_prefixes:
            return None
        else:
            return attributes[self.attrib]


class ToplevelOrthoXMLParser(object):
    def __init__(self, gene_mapper):
        self.gene_mapper = gene_mapper
        self.cur_hog_id = None
        self.cur_hog_depth = 0
        self.cur_hog_memb = []
        self.genes = {}
        self.levels = {}
        self.hogs = {}

    def start(self, tag, attrib):
        if tag == "{http://orthoXML.org/2011/}gene":
            self.genes[int(attrib["id"])] = self.gene_mapper.map(attrib)
        elif tag == "{http://orthoXML.org/2011/}geneRef":
            gene_id = int(attrib["id"])
            gene = self.genes[gene_id]
            if gene is not None:
                self.cur_hog_memb.append(gene)
        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            if self.cur_hog_depth == 0:
                self.cur_hog_id = int(attrib["id"])
            self.cur_hog_depth += 1
        elif (
            tag == "{http://orthoXML.org/2011/}property"
            and attrib["name"] == "TaxRange"
        ):
            if self.cur_hog_depth == 1:
                self.levels[self.cur_hog_id] = attrib["value"]

    def end(self, tag):
        if tag == "{http://orthoXML.org/2011/}orthologGroup":
            self.cur_hog_depth -= 1
            if self.cur_hog_depth == 0:
                self.hogs[self.cur_hog_id] = frozenset(self.cur_hog_memb)
                self.cur_hog_memb = []

    def data(self, data):
        pass

    def close(self):
        return


class Analysis(object):
    def __init__(self, hogs, rootlevels):
        self.hogs = hogs
        self.rootlevels = rootlevels


def load_toplevel_hogs(fpath, species_to_skip=None):
    with auto_open(fpath) as fh:
        if species_to_skip is None:
            gene_mapper = GeneMapper()
        else:
            gene_mapper = FilteredSpeciesGeneMapper(species_to_skip)
        toplevelparser = ToplevelOrthoXMLParser(gene_mapper)
        parser = XMLParser(target=toplevelparser)
        for line in fh:
            parser.feed(line)
    return Analysis(toplevelparser.hogs, toplevelparser.levels)
