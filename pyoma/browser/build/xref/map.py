import abc
import collections
import os
import pickle
from typing import Mapping, Set, Union, List
import logging

import networkx as nx
import tables
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


from ...tablefmt import XRefTable
from ...db import SequenceSearch
from ...db import Database, OmaIdMapper
from ....common import auto_open
from ..builder import XrefStorer

logger = logging.getLogger(__name__)
TaxRange = collections.namedtuple("TaxRange", ["genomes", "entry_nr_range"])
Match = collections.namedtuple("Match", ["id", "entries", "method", "identity"])


class Mapper(metaclass=abc.ABCMeta):
    """Mapper of Bio.SeqRecords to proteins based on sequence"""

    identity_threshold = 0.9

    def __init__(
        self,
        db_path: os.PathLike,
        seq_idx_path: os.PathLike,
        src_xref_path: os.PathLike,
        taxid_mapping: Mapping[int, Set[int]],
    ):
        self.db = Database(db_path)
        self.searcher = SequenceSearch(self.db, seq_idx_path)
        self.oma_id_mapper: OmaIdMapper = self.db.id_mapper["OMA"]
        self.taxid_mapping = self._identify_entry_ranges_for_taxid_mappings(taxid_mapping)
        self.src_xrefs = self._load_source_ids(src_xref_path)

    @abc.abstractmethod
    def get_taxid(self, rec):
        pass

    def _identify_entry_ranges_for_taxid_mappings(
        self, taxid_mappings: Mapping[int, Set[int]]
    ) -> Mapping[int, TaxRange]:
        res = {}
        for src_taxid, target_taxids in taxid_mappings.items():
            ranges, genomes = [], []
            for taxid in target_taxids:
                g = self.db.tax.genomes[taxid]
                ranges.append((g.entry_nr_offset + 1, g.entry_nr_offset + g.nr_entries))
                genomes.append(g)
            ranges = sorted(ranges, key=lambda x: x[0])
            for k in range(len(ranges) - 1):
                assert (
                    ranges[k][1] + 1 == ranges[k + 1][0]
                ), f"ranges not as expected for {src_taxid} -> {target_taxids}: {ranges}"
            res[src_taxid] = TaxRange(genomes, (ranges[0][0], ranges[-1][1]))
        return res

    def _load_source_ids(self, dbpath: os.PathLike) -> Mapping[str, Set[int]]:
        mapping = collections.defaultdict(set)
        with tables.open_file(dbpath) as xf:
            for row in xf.root.XRef:
                mapping[row["XRefId"].decode()].add(int(row["EntryNr"]))
        return mapping

    def check_xrefs(self, rec):
        res = set()
        avoid_prefix = {
            "GO",
            "OMA",
            "KEGG",
            "Araport",
            "eggNOG",
            "InParanoid",
            "OrthoDB",
            "PhylomeDB",
            "PRO",
            "Proteomes",
            "ExpressionAtlas",
            "Gene3D",
            "InterPro",
            "PANTHER",
            "Pfam",
            "PIRSF",
            "PRINTS",
            "SMART",
            "SUPFAM",
            "PROSITE",
            "BioGRID",
            "IntAct",
            "AlphaFoldDB",
            "SMR",
        }
        if rec.id in self.src_xrefs:
            res.update(self.src_xrefs[rec.id])
        for xref in rec.dbxrefs:
            prefix, id_ = xref.split(":", maxsplit=1)
            if prefix not in avoid_prefix:
                if id_ in self.src_xrefs:
                    res.update(self.src_xrefs[id_])
        return res

    def map_record(self, rec: SeqRecord) -> Union[None, Match]:
        """Map a SeqRecord to the database by comparing the sequence with the proteins in OMA.

        The following steps are done:

          - minimum length of sequence (10 AA)
          - annotated taxid is relevant for OMA dataset
          - exact match of sequence (within relevant genomes)
          - if exact match found, return these matches. In case >1, we first check if
            one of the db_xrefs has one of the OMA source IDs/ACs listed.
          - if no exact match found, do an approximate match search, If the most similar
            one (highest score in global alignment) exceeds the identity threshold, this
            match is returned.
        """
        # minimal sequence length
        if len(rec.seq) < 10:
            logger.info(f"skipping {rec.id} (too short: {len(rec.seq)} AA)")
            return None
        # check that taxid of xref is among relevant taxids
        xref_taxid = self.get_taxid(rec)
        if xref_taxid not in self.taxid_mapping:
            return None

        # exact match, allowing for 1 AA difference in seq length (often added start codon)
        taxrange = self.taxid_mapping[xref_taxid]
        entry_nrs = set(
            self.searcher.exact_search(
                str(rec.seq), only_full_length=True, max_len_diff=1, entrynr_range=taxrange.entry_nr_range
            )
        )
        logger.debug(f"{rec.id} maps to {len(entry_nrs)} exactly")
        if len(entry_nrs) > 1 and len(taxrange.genomes) == 1:
            src_xref_match = self.check_xrefs(rec)
            seq_and_id_support = src_xref_match.intersection(entry_nrs)
            logger.debug(
                f"{rec.id} contains {len(src_xref_match)} crossreferences, of which {len(seq_and_id_support)} are in the seq set."
            )
            if seq_and_id_support:
                return Match(rec.id, seq_and_id_support, "exact", 1)
        if len(entry_nrs) > 0:
            return Match(rec.id, entry_nrs, "exact", 1)

        # now, we try approximate search
        approx_matches = self.searcher.approx_search(
            str(rec.seq), n=20, entrynr_range=taxrange.entry_nr_range, coverage=0.7, alignment="global"
        )
        logger.debug(f"computed global alignments for {len(approx_matches)} approx matches")
        if len(approx_matches) > 0:
            s1, s2 = (approx_matches[0][1]["alignment"][0][0], approx_matches[0][1]["alignment"][1][0])
            identity = sum(1 for b1, b2 in zip(s1, s2) if b1 == b2) / len(s1)
            logger.info(
                f"best approximate match of {rec.id} is entry_nr {approx_matches[0][0]}: "
                f"{identity:.3f} identy, score {approx_matches[0][1]['score']:.3f}"
            )
            logger.debug(f"Alignment:\n{s1}\n{s2}")
            if identity > self.identity_threshold:
                return Match(rec.id, {approx_matches[0][0]}, "approx", identity)
        return None


class SwissFormatMapper(Mapper):
    def get_taxid(self, rec):
        return int(rec.annotations.get("ncbi_taxid")[0])


class SwissProtMapper(SwissFormatMapper):
    identity_threshold = 0.6


class GenbankFormatMapper(Mapper):
    def get_taxid(self, rec):
        try:
            db_xrefs = rec.features[0].qualifiers["db_xref"]
            for ref in db_xrefs:
                key, val = ref.split(":", 1)
                if key == "taxon":
                    return int(val)
        except (KeyError, IndexError, AttributeError, ValueError):
            logger.warning("cannot extract taxid from %s", rec)
        return 0


def map_xrefs(
    fpath: os.PathLike,
    format: str,
    source: str,
    out_fpath: os.PathLike,
    db: os.PathLike,
    seq_idx: os.PathLike,
    xref_db: os.PathLike,
    taxid_mapping: Mapping[int, Set[int]],
):
    mapper_cls = (
        SwissProtMapper if source == "swissprot" else SwissFormatMapper if source == "trembl" else GenbankFormatMapper
    )
    mapper = mapper_cls(db, seq_idx, xref_db, taxid_mapping)
    mapping_results = []
    with auto_open(fpath, "rt") as fh:
        rec_iter = SeqIO.parse(fh, format=format)
        for rec in rec_iter:
            res = mapper.map_record(rec)
            if res is not None:
                mapping_results.append(res)
    with auto_open(out_fpath, "wb") as fh:
        pickle.dump(mapping_results, fh)
    logger.info(f"wrote {len(mapping_results)} records to {out_fpath}")


def _filter_graph(G):
    def graph_to_tuples(G):
        for u, v, data in G.edges(data=True):
            if isinstance(u, int):
                u, v = v, u
            yield u, v, data["weight"], data.get("propagate", True)

    # Numeric nodes ==> entry nr in OMA; string nodes ==> source ids
    for cc in nx.connected_components(G):
        SG = G.subgraph(cc).copy()
        if len(cc) <= 2:
            yield from graph_to_tuples(SG)
        else:
            src = {n for n in SG.nodes() if isinstance(n, str)}
            tar = {n for n in SG.nodes() if isinstance(n, int)}
            if len(src) == 1 or len(tar) == 1:
                # single source node maps to several sequences. We keep all of them
                # (should all have the same maximal weight)
                # several source ids map to only one oma entry. Keep only the top-matching
                # # with propagation and for the other only the
                max_sim = None
                for u, v, data in sorted(SG.edges(data=True), key=lambda e: -e[2]["weight"]):
                    if max_sim is None:
                        max_sim = data["weight"]
                    if data["weight"] < max_sim * 0.8:
                        data["propagate"] = False
                yield from graph_to_tuples(SG)
            else:
                # we have a NxM case. let's compute a maximal matching and keep only those links
                max_matching = nx.max_weight_matching(SG, maxcardinality=True, weight="weight")
                for x in max_matching:
                    M = SG.subgraph(x).copy()
                    yield from graph_to_tuples(M)


def identify_best_matching(map_files: List[os.PathLike], out: os.PathLike):
    G = nx.Graph()
    for map_file in map_files:
        with open(map_file, "rb") as fh:
            matches = pickle.load(fh)
        for match in matches:
            edges = [(match.id, z, 2 if match.method == "exact" else match.identity) for z in match.entries]
            G.add_weighted_edges_from(edges)

    final_matches = [_filter_graph(G)]
    with auto_open(out, "wb") as fh:
        pickle.dump(final_matches, fh)
