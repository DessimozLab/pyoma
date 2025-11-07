import collections
import io
import logging
import re
from pathlib import Path
from typing import Set, List, Mapping, Union
import ete3
import omataxonomy

from ....common import auto_open

logger = logging.getLogger(__name__)


class ChunkWriter:
    def __init__(self, prefix, chunksize):
        self.prefix = Path(prefix)
        self.chunksize = chunksize
        self.chunk = 0
        self.in_chunk = 0
        self._fp = None

    def __enter__(self):
        self._open_new_chunk()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._fp:
            self._fp.close()

    def _open_new_chunk(self):
        self.chunk += 1
        fname = f"{self.prefix.stem}-{self.chunk:04d}.gz"
        fpath = self.prefix.parent / fname
        if self._fp:
            self._fp.close()
        logger.info(f"Opening new chunk: {fpath}")
        self._fp = auto_open(fpath, "wb")

    def write(self, chunk: bytes):
        self._fp.write(chunk)
        self.in_chunk += 1
        if self.in_chunk >= self.chunksize:
            self._open_new_chunk()
            self.in_chunk = 0


class SwissSpeciesFilter(object):
    RECORD_SEP = b"//\n"
    TAXON_RE = re.compile(rb"^OX\s+NCBI_TaxID=(?P<taxid>\d+)")
    TAXON_LINE_STR = b"NCBI_TaxID="

    def __init__(self, fpath, taxids):
        self.fpath = fpath
        self.tax_of_interest = frozenset(taxids)

    def __enter__(self):
        self.file = auto_open(self.fpath, "rb")
        self.buffer = io.BufferedReader(self.file)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.buffer:
            self.buffer.close()
        if self.file:
            self.file.close()

    def filter(self, chunk_writer, log_every=10_000):
        """yields the sequence records of taxons of interest"""
        record_chunks = []
        taxid = None
        record_count = 0
        for line in self.buffer:
            record_chunks.append(line)
            if self.TAXON_LINE_STR in line:
                m = self.TAXON_RE.match(line)
                if m:
                    taxid = int(m.group("taxid"))
                    if taxid not in self.tax_of_interest:
                        taxid = 0
            if line == self.RECORD_SEP:
                if taxid in self.tax_of_interest:
                    chunk_writer.write(b"".join(record_chunks))
                record_chunks.clear()
                taxid = None

                record_count += 1
                if record_count % log_every == 0:
                    logger.info("Processed %d records from %s", record_count, self.fpath)


class GenbankSpeciesFilter(SwissSpeciesFilter):
    TAXON_RE = re.compile(rb'^\s+/db_xref="taxon:(?P<taxid>\d+)"')
    TAXON_LINE_STR = b'/db_xref="taxon:'


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
    logger.debug(f"mapped {len(species_taxids)} taxids to {len(ranks)} ranks: {ranks}")
    relevant_taxids = collections.defaultdict(set)
    for taxid in species_taxids:
        relevant_taxids[taxid].add(taxid)
    # for every species (that has a limited rank), select all its children (recursively)
    for taxid in species_taxids:
        if ranks.get(taxid) in ("species", "genus", "varietas", "strain"):
            sub_taxids = ncbi_taxonomy.get_descendant_taxa(taxid, intermediate_nodes=True)
            for sub_taxid in sub_taxids:
                relevant_taxids[sub_taxid].add(taxid)
    logger.debug(f"relevant taxid mapping: {relevant_taxids}")
    # now, build a tree with the input species and select the genus rank nodes.
    # for each of those, store in the mapping every subnode taxid to a all the
    # species nodes in that clade.
    ncbi_taxids = list(ncbi_taxonomy.get_taxid_translator(species_taxids).keys())
    if len(ncbi_taxids) == 0:
        return relevant_taxids

    tree = ncbi_taxonomy.get_topology(ncbi_taxids, intermediate_nodes=True, annotate=True)
    for genus_node in tree.iter_search_nodes(rank="genus"):
        for nn in genus_node.traverse(strategy="postorder"):
            if nn.is_leaf():
                nn.add_feature("subtaxids", {nn.taxid})
            else:
                nn.add_feature("subtaxids", set.union(*list(x.subtaxids for x in nn.get_children())))
                relevant_taxids[nn.taxid].update(nn.subtaxids)
    logger.debug(f"final relevant taxid mapping: {relevant_taxids}")
    return relevant_taxids
