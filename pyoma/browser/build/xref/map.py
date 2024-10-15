import collections
import os
from typing import Mapping, Set, Tuple
import logging

import numpy
import tables
from Bio import SeqIO

from ...tablefmt import XRefTable
from ...db import SequenceSearch
from ...db import Database, OmaIdMapper
from ....common import auto_open
from ..builder import XrefStorer

logger = logging.getLogger(__name__)
TaxRange = collections.namedtuple("TaxRange", ["genomes", "entry_nr_range"])


class Mapper:
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

    def get_taxid(self, rec):
        return int(rec.annotations.get("ncbi_taxid")[0])

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
                assert ranges[k][1] + 1 == ranges[k + 1][0], f"ranges not as expected for {src_taxid} -> {target_taxids}: {ranges}"
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
        for xref in rec.dbxrefs:
            prefix, id_ = xref.split(":", maxsplit=1)
            if prefix not in avoid_prefix:
                if id_ in self.src_xrefs:
                    res.update(self.src_xrefs[id_])
        return res

    def map_record(self, rec):
        # minimal sequence length
        if len(rec.seq) < 10:
            return
        # check that taxid of xref is among relevant taxids
        xref_taxid = self.get_taxid(rec)
        if xref_taxid not in self.taxid_mapping:
            return

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
                return seq_and_id_support, "exact"
        if len(entry_nrs) > 0:
            return entry_nrs, "exact"

        # now, we try approximate search
        approx_matches = self.searcher.approx_search(
            str(rec.seq), n=20, entrynr_range=taxrange.entry_nr_range, coverage=0.7, alignment="global"
        )
        logger.debug(f"computed global alignments for {len(approx_matches)} approx matches")
        if len(approx_matches) > 0:
            s1, s2 = (approx_matches[0][1]["alignment"][0][0], approx_matches[0][1]["alignment"][1][0])
            identity = sum(1 for b1, b2 in zip(s1, s2) if b1 == b2) / len(s1)
            logger.info(
                f"best approximate match of {rec.id} is entry_nr {approx_matches[0][0]}: {identity:.3f} identy, score {approx_matches[0][1]['score']:.3f}"
            )
            logger.debug(f"Alignment:\n{s1}\n{s2}")
            if identity > 0.90:
                return {approx_matches[0][0]}, "approx"
        return set([])


def map_xrefs(
    fpath: os.PathLike,
    format: str,
    db: os.PathLike,
    seq_idx: os.PathLike,
    xref_db: os.PathLike,
    taxid_mapping: Mapping[int, Set[int]],
):
    mapper = Mapper(db, seq_idx, xref_db, taxid_mapping)
    with auto_open(fpath, "rt") as fh:
        rec_iter = SeqIO.parse(fh, format=format)
        with open("xref.h5", "wt") as fout:
            for rec in rec_iter:
                res = mapper.map_record(rec)
                fout.write(str(res))
                fout.write("\n")
