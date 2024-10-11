import collections
import logging
import re
from typing import Set, List, Mapping, Union
import ete3
import omataxonomy

from ....common import auto_open

logger = logging.getLogger(__name__)


class ChunkWriter:
    def __init__(self, prefix, chunksize):
        self.prefix = prefix
        self.chunksize = chunksize

    def __enter__(self):
        self.chunk = 1
        self.in_chunk = 0
        self._fp = auto_open(f"{self.prefix}-{self.chunk:03d}.gz", "wb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._fp.close()

    def write(self, chunk: bytes):
        self._fp.write(chunk)
        self.in_chunk += 1
        if self.in_chunk >= self.chunksize:
            self._fp.close()
            self.chunk += 1
            self.in_chunk = 0
            self._fp = auto_open(f"{self.prefix}-{self.chunk:03d}.gz", "wb")


class SwissSpeciesFilter(object):
    RECORD_SEP = b"//\n"
    TAXON_RE = re.compile(rb"^OX\s+NCBI_TaxID=(?P<taxid>\d+)")

    def __init__(self, fpath, taxids):
        self.fpath = fpath
        self.tax_of_interest = frozenset(taxids)

    def __enter__(self):
        self.file = auto_open(self.fpath, "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            self.file.close()

    def filter(self, chunk_writer):
        """yields the sequence records of taxons of interest"""
        self.file.seek(0)
        os, buf = 0, []
        for line_nr, line in enumerate(self.file):
            buf.append(line)
            m = self.TAXON_RE.match(line)
            if m:
                os = int(m.group("taxid"))
                if os not in self.tax_of_interest:
                    os = 0
            elif line == self.RECORD_SEP:
                if os != 0:
                    chunk_writer.write(b"".join(buf))
                buf = []
            if line_nr % 1e6 == 0:
                logger.info("processed %d lines from %s", line_nr, self.fpath)


class GenbankSpeciesFilter(SwissSpeciesFilter):
    TAXON_RE = re.compile(rb'^\s+/db_xref="taxon:(?P<taxid>\d+)"')


def filter_records_on_taxids(fpath: str, writer, taxids: Set[int], format: str):
    if format == "swiss":
        handler = SwissSpeciesFilter
    elif format == "genbank":
        handler = GenbankSpeciesFilter
    else:
        raise ValueError("format must be either 'swiss' or 'genbank'")

    with handler(fpath, taxids) as parser:
        parser.filter(chunk_writer=writer)


def load_relevant_taxids(
    species_taxids, ncbi_taxonomy: Union[ete3.NCBITaxa, omataxonomy.Taxonomy]
) -> Mapping[int, Set[int]]:
    """load relevant taxids for xref mapping for a given set of species.

    relevant taxids are the sub-taxids of the species selected, and
    also the parent taxids up to the genus rank (Higher parents are
    left out, as they are likely too general)"""
    ranks = ncbi_taxonomy.get_rank(species_taxids)
    relevant_taxids = collections.defaultdict(set)
    for taxid in species_taxids:
        relevant_taxids[taxid].add(taxid)
    # for every species (that has a limited rank), select all its children (recursively)
    for taxid in species_taxids:
        if ranks[taxid] in ("species", "genus", "varietas", "strain"):
            sub_taxids = ncbi_taxonomy.get_descendant_taxa(taxid, intermediate_nodes=True)
            for sub_taxid in sub_taxids:
                relevant_taxids[sub_taxid].add(taxid)
    # now, build a tree with the input species and select the genus rank nodes.
    # for each of those, store in the mapping every subnode taxid to a all the
    # species nodes in that clade.
    tree = ncbi_taxonomy.get_topology(species_taxids, intermediate_nodes=True, annotate=True)
    for genus_node in tree.iter_search_nodes(rank="genus"):
        for nn in genus_node.traverse(strategy="postorder"):
            if nn.is_leaf():
                nn.add_feature("subtaxids", {nn.taxid})
            else:
                nn.add_feature("subtaxids", set.union(*list(x.subtaxids for x in nn.get_children())))
                relevant_taxids[nn.taxid].update(nn.subtaxids)
    return relevant_taxids
