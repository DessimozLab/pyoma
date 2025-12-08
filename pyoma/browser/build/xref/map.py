import abc
import collections
import heapq
import itertools
import os
import pickle
import re
from concurrent.futures import ProcessPoolExecutor, as_completed

from typing import Mapping, Set, Union, List, Tuple, Callable
import logging

import networkx as nx
import pandas as pd
import tables
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from ...db import SequenceSearch
from ...db import Database, OmaIdMapper
from ....common import auto_open
from ..builder import XrefStorer

logger = logging.getLogger(__name__)
TaxRange = collections.namedtuple("TaxRange", ["genomes", "entry_nr_range"])
Match = collections.namedtuple("Match", ["id", "entries", "method", "identity"])
BestMatch = collections.namedtuple("BestMatch", ["id", "entry", "identity", "propagate"])


class Mapper(metaclass=abc.ABCMeta):
    """Mapper of Bio.SeqRecords to proteins based on sequence"""

    identity_threshold = 0.9

    def __init__(
        self,
        db_path: os.PathLike,
        seq_idx_path: os.PathLike,
        src_xref_path: os.PathLike,
        taxid_mapping: Mapping[int, Set[int]],
        approx_align: bool = True,
    ):
        self.db = Database(db_path)
        self.searcher = SequenceSearch(self.db, seq_idx_path)
        self.oma_id_mapper: OmaIdMapper = self.db.id_mapper["OMA"]
        self.taxid_mapping = self._identify_entry_ranges_for_taxid_mappings(taxid_mapping)
        self.src_xrefs = self._load_source_ids(src_xref_path)
        self._do_approx_align = approx_align

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
            if all(ranges[k][1] + 1 == ranges[k + 1][0] for k in range(len(ranges) - 1)):
                tax_range = TaxRange(genomes, (ranges[0][0], ranges[-1][1]))
            else:
                tax_range = TaxRange(genomes, set(r for rng in ranges for r in range(rng[0], rng[1] + 1)))
            res[src_taxid] = tax_range
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

        if self._do_approx_align:
            # now, we try approximate search
            approx_matches = self.searcher.approx_search(
                str(rec.seq),
                n=20,
                entrynr_range=taxrange.entry_nr_range,
                coverage=self.identity_threshold - 0.2,
                alignment="global",
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
        else:
            kmer_matches = self.searcher.approx_search_no_align(
                str(rec.seq), coverage=self.identity_threshold - 0.2, entrynr_range=taxrange.entry_nr_range
            )
            logger.debug(f"searched for kmer-based matches: {len(kmer_matches)} approx matches")
            if len(kmer_matches) > 1:
                logger.debug(
                    f"  -> {kmer_matches[0][1]} vs {kmer_matches[1][1]}: {kmer_matches[0][1]/kmer_matches[1][1]:.3f} ratio best/second"
                )
            return Match(rec.id, {kmer_matches[0][0]}, "approx", kmer_matches[0][1]) if kmer_matches else None
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


def map_chunk_of_xrefs_worker(args):
    mapper_cls, db, seq_idx, xref_db, taxid_mapping, align, recs = args
    mapper = mapper_cls(db, seq_idx, xref_db, taxid_mapping, align)
    mapping_results = []
    for rec in recs:
        res = mapper.map_record(rec)
        if res is not None:
            mapping_results.append(res)
    return mapping_results


def map_xrefs(
    fpaths: List[os.PathLike],
    format: str,
    source: str,
    out_fpath: os.PathLike,
    db: os.PathLike,
    seq_idx: os.PathLike,
    xref_db: os.PathLike,
    taxid_mapping: Mapping[int, Set[int]],
    align: bool = True,
    nr_procs: int = 1,
):
    mapper_cls = (
        SwissProtMapper if source == "swissprot" else SwissFormatMapper if source == "trembl" else GenbankFormatMapper
    )
    logger.info(f"mapping {len(fpaths)} xref files using {mapper_cls.__name__} with {nr_procs} processes")

    def chunkify(fpaths, size=50):
        chunk = []
        for fpath in fpaths:
            with auto_open(fpath, "rt") as fh:
                rec_iter = SeqIO.parse(fh, format=format)
                for rec in rec_iter:
                    chunk.append(rec)
                    if len(chunk) == size:
                        yield chunk
                        chunk = []
        if chunk:
            yield chunk

    with ProcessPoolExecutor(max_workers=nr_procs) as pool:
        futures = []
        for chunk in chunkify(fpaths, size=50):
            args = (mapper_cls, db, seq_idx, xref_db, taxid_mapping, align, chunk)
            futures.append(pool.submit(map_chunk_of_xrefs_worker, args))

        mapping_results = []
        for fut in as_completed(futures):
            res = fut.result()
            if res:
                mapping_results.extend(res)
    with open(out_fpath, "wb") as fh:
        pickle.dump(mapping_results, fh)
    logger.info(f"wrote {len(mapping_results)} records to {out_fpath}")


def _filter_graph(G):
    def graph_to_tuples(G):
        for u, v, data in G.edges(data=True):
            if isinstance(u, int):
                u, v = v, u
            yield BestMatch(u, v, data["weight"], data.get("propagate", True))

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


def identify_best_matching(map_files: List[os.PathLike]) -> Mapping[str, List[BestMatch]]:
    G = nx.Graph()
    for map_file in map_files:
        with open(map_file, "rb") as fh:
            matches = pickle.load(fh)
        for match in matches:
            edges = [(match.id, z, 2 if match.method == "exact" else match.identity) for z in match.entries]
            G.add_weighted_edges_from(edges)

    final_matches = collections.defaultdict(list)
    for best_match in _filter_graph(G):
        final_matches[best_match.id].append(best_match)
    return final_matches


class CrossRefsExtractor(metaclass=abc.ABCMeta):
    SRC_ENUM_KEY = "SourceID"

    def __init__(self, out_fpath: os.PathLike, match_lookup: Mapping[str, List[BestMatch]]):
        self.storer = XrefStorer(out_fpath, buffer_size=10_000_000, multi_files=True)
        self.match_lookup = match_lookup
        self.src_enum_val = None

    def __enter__(self):
        self.storer.__enter__()
        self.src_enum_val = self.storer.source_enum[self.SRC_ENUM_KEY]
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.storer.__exit__(exc_type, exc_val, exc_tb)

    @abc.abstractmethod
    def extract_crossrefs(self, rec: SeqRecord) -> Tuple[List[Tuple[int, str]], List[Tuple[int, str]], List[str]]:
        pass

    def map_record(self, rec: SeqRecord) -> None:
        def score2verif(score):
            if score < 1:
                return self.storer.verify_enum["modified"]
            return self.storer.verify_enum["exact"]

        if rec.id in self.match_lookup:
            crossrefs, propagated_crossrefs, ec = self.extract_crossrefs(rec)
            for match in self.match_lookup[rec.id]:
                verif = score2verif(match.identity)
                ident = min(match.identity, 1)
                for src, xrefid in crossrefs:
                    self.storer.add_xref(match.entry, src, xrefid, verif, ident)
                if match.propagate:
                    verif_prop = max(verif, self.storer.verify_enum["unchecked"])
                    for src, xrefid in propagated_crossrefs:
                        self.storer.add_xref(match.entry, src, xrefid, verif_prop, ident)
                    for ec_term in ec:
                        self.storer.add_ec(match.entry, ec_term)


class UniProtKBCrossRefsExtractor(CrossRefsExtractor):
    SRC_ENUM_KEY = "UniProtKB/TrEMBL"
    PROT_NAME_RE = re.compile(r"^(?P<typ>((Rec)|(Alt)|(Sub))Name): Full=(?P<name>[^{]*)")
    ENS_RE = re.compile(r"ENS(?P<species>[A-Z]{0,3})(?P<typ>[GTP])(?P<num>\d{11})")

    def __enter__(self):
        super().__enter__()
        key_map = {
            "Name": "Gene Name",
            "Synonyms": "Synonym",
            "OrderedLocusNames": "Ordered Locus Name",
            "ORFNames": "ORF Name",
            "HGNC": "HGNC",
            "RecName": "Protein Name",
            "SubName": "Protein Name",
            "AltName": "Alternative Protein Name",
            "EnsemblPlants": "EnsemblGenomes",
            "EnsemblFungi": "EnsemblGenomes",
            "EnsemblMetazoa": "EnsemblGenomes",
            "EnsemblProtists": "EnsemblGenomes",
            "EnsemblBacteria": "EnsemblGenomes",
            "Bgee": "Bgee",
            "SMR": "Swiss Model",
            "PDB": "PDB",
            "STRING": "STRING",
            "neXtProt": "neXtProt",
            "EPD": "EPD",
            "GlyConnect": "GlyConnect",
            "GeneID": "EntrezGene",
            "WikiGene": "WikiGene",
            "RefSeq": "RefSeq",
            "KEGG": "KEGG",
            "AGR": "AGR",
        }
        self.key_map = {k: self.storer.source_enum[v] for k, v in key_map.items()}
        return self

    def iter_gene_names(self, annotations: Mapping) -> List[Tuple[int, str]]:
        """extract the gene names (Name, Synonyms, OrderedLocusNames, ORFNames)"""
        try:
            gene_names = annotations["gene_name"]
        except KeyError:
            return
        for elem in gene_names:
            for typ, value in elem.items():
                try:
                    typ = self.key_map[typ]
                except KeyError:
                    continue
                value = [value] if isinstance(value, str) else value
                for val in value:
                    yield typ, val

    def iter_parsed_description(self, desc: str):
        chunks = desc.split("; ")
        for chunk in chunks:
            if m := self.PROT_NAME_RE.match(chunk):
                typ = self.key_map[m.group("typ")]
                val = m.group("name").strip()
                yield typ, val
            elif chunk.startswith("EC="):
                ec = chunk[3:].split(" ")[0]
                yield "ec", ec

    def iter_crossrefs(self, xrefs):
        for xref in xrefs:
            db, ref = xref.split(":", maxsplit=1)
            try:
                typ = self.key_map[db]
                yield typ, ref
            except KeyError:
                if db == "Ensembl":
                    m = self.ENS_RE.match(ref)
                    if m:
                        if m.group("typ") == "P":
                            typ = self.storer.source_enum["Ensembl Protein"]
                        elif m.group("typ") == "G":
                            typ = self.storer.source_enum["Ensembl Gene"]
                        else:
                            typ = self.storer.source_enum["Ensembl Transcript"]
                        yield typ, ref

    def extract_crossrefs(self, rec):
        stable, projected, ec = [], [], []
        stable.append((self.src_enum_val, rec.id))
        try:
            stable.append((self.src_enum_val, rec.name))
        except AttributeError:
            pass
        # extract gene names
        stable.extend(self.iter_gene_names(rec.annotations))
        # extract protein names and ec
        for typ, val in self.iter_parsed_description(rec.description):
            if typ == "ec":
                ec.append(val)
            else:
                stable.append((typ, val))
        # add crossreferences (as projected)
        projected.extend(self.iter_crossrefs(rec.dbxrefs))
        return stable, projected, ec


class SwissProtCrossRefsExtractor(UniProtKBCrossRefsExtractor):
    SRC_ENUM_KEY = "UniProtKB/SwissProt"


class RefSeqCrossRefsExtractor(CrossRefsExtractor):
    SRC_ENUM_KEY = "RefSeq"
    GENEID_RE = re.compile(r"^GeneID:(?P<id>\d+)$")

    def __enter__(self):
        super().__enter__()
        self.geneid_val = self.storer.source_enum["EntrezGene"]
        return self

    def extract_crossrefs(self, rec):
        crossrefs = []
        crossrefs.append((self.src_enum_val, rec.id))
        for f in rec.features:
            if f.type == "CDS":
                try:
                    gene_id = f.qualifiers["db_xref"][0]
                    m = self.GENEID_RE.match(gene_id)
                    if m:
                        crossrefs.append((self.geneid_val, m.group("id")))
                except (KeyError, AttributeError, IndexError):
                    pass
        return crossrefs, [], []


def collect_crossrefs(
    xrefs: List[os.PathLike], source: str, format: str, map_files: List[os.PathLike], out: os.PathLike
):
    best_matches = identify_best_matching(map_files)
    if source == "swissprot":
        collector_cls = SwissProtCrossRefsExtractor
    elif source == "trembl":
        collector_cls = UniProtKBCrossRefsExtractor
    elif source == "refseq":
        collector_cls = RefSeqCrossRefsExtractor
    else:
        raise ValueError(f"Unknown source: {source}")
    with collector_cls(out, best_matches) as collector:
        for fpath in xrefs:
            with auto_open(fpath, "rt") as fh:
                rec_iter = SeqIO.parse(fh, format)
                for rec in rec_iter:
                    collector.map_record(rec)


def merge_sorted_h5_table(
    input_h5_handles: List[tables.File],
    table_path: str,
    sort_columns: List[str],
    dupl_subset_columns: List[str],
    storer_callback: Callable,
    enr_batch_size: int = 100,
    chunksize: int = 500_000,
):
    """
    Merge pre-sorted HDF5 tables in memory-bounded, batch-wise fashion,
    deduplicating per EntryNr across files.
    """

    # Build iterators for each file
    def iter_sorted_with_carryover(h5: tables.File):
        carryover = pd.DataFrame()
        table = h5.get_node(table_path)
        nrows = table.nrows
        colnames = table.colnames
        for start in range(0, nrows, chunksize):
            stop = min(start + chunksize, nrows)
            chunk = pd.DataFrame.from_records(table.read(start, stop), columns=colnames)
            if not carryover.empty:
                chunk = pd.concat([carryover, chunk], ignore_index=True)
                carryover = pd.DataFrame()
            if not chunk.empty:
                last_enr = chunk["EntryNr"].iloc[-1]
                mask = chunk["EntryNr"] == last_enr
                carryover = chunk.loc[mask]
                chunk = chunk.loc[~mask]
            if not chunk.empty:
                yield chunk
        if not carryover.empty:
            yield carryover

    # store dtype for later use
    dt = input_h5_handles[0].get_node(table_path).dtype
    dt = {c: dt[c] for c in dt.names}

    readers = [iter_sorted_with_carryover(path) for path in input_h5_handles]
    current_chunks = [next(r, pd.DataFrame()) for r in readers]

    while any(not df.empty for df in current_chunks):
        # Find the smallest available EntryNr across all files
        min_enr = min(df["EntryNr"].iloc[0] for df in current_chunks if not df.empty)
        enr_limit = min_enr + enr_batch_size

        per_batch = []
        for i, df in enumerate(current_chunks):
            if df.empty:
                continue
            # Take all rows < enr_limit
            mask = df["EntryNr"] < enr_limit
            batch_df = df.loc[mask]
            leftover_df = df.loc[~mask]
            if not batch_df.empty:
                per_batch.append(batch_df)
                # Load next chunk if necessary
            current_chunks[i] = leftover_df if not leftover_df.empty else next(readers[i], pd.DataFrame())

        if not per_batch:
            continue

        # Combine batch and find the true highest EntryNr seen
        merged = pd.concat(per_batch, ignore_index=True)
        max_enr_in_batch = merged["EntryNr"].max()

        # Add any leftover rows from files that belong to the last EntryNr
        extra_rows = []
        for i, df in enumerate(current_chunks):
            if not df.empty and df["EntryNr"].iloc[0] == max_enr_in_batch:
                mask = df["EntryNr"] == max_enr_in_batch
                extra_rows.append(df.loc[mask])
                current_chunks[i] = df.loc[~mask]

        if extra_rows:
            merged = pd.concat([merged, *extra_rows], ignore_index=True)

        merged.sort_values(by=sort_columns, inplace=True)
        merged.drop_duplicates(subset=dupl_subset_columns, keep="first", inplace=True)

        # Pass final merged batch to storer callback
        storer_callback(merged.to_records(index=False))
        logger.debug("Processed up to EntryNr %s", max_enr_in_batch)
    logger.info("Finished merging table %s", table_path)


def _fetch_combine_and_reduce_input_data(h5_handles, table_path, sort_columns, dupl_subset_columns, storer_callback):
    dt = h5_handles[0].get_node(table_path).dtype
    dt = {c: dt[c] for c in dt.names}
    iters = [
        map(lambda row: row.fetch_all_fields(), h5.get_node(table_path).itersorted(sortby="EntryNr"))
        for h5 in h5_handles
    ]
    a_tab = h5_handles[0].get_node(table_path)
    try:
        tot_entries = a_tab[a_tab.colindexes["EntryNr"][-1]]["EntryNr"]
    except IndexError:
        tot_entries = len(a_tab)  # most likely empty table
    queue = heapq.merge(*iters, key=lambda row: row["EntryNr"])
    for enr, data_per_enr_it in tqdm(itertools.groupby(queue, key=lambda row: row["EntryNr"]), total=tot_entries):
        df = pd.DataFrame.from_records(data_per_enr_it)
        df.sort_values(by=sort_columns, inplace=True)
        df.drop_duplicates(subset=dupl_subset_columns, keep="first", inplace=True)
        storer_callback(df.to_records(index=False, column_dtypes=dt).tolist())


def combine_xrefs(xrefs: List[os.PathLike], out: os.PathLike):
    with XrefStorer(
        str(out), index_cols=["EntryNr", "XRefId", "XRefSource"], suffix_col="XRefId", multi_files=False
    ) as storer:
        h5hs = [tables.open_file(str(fn), mode="r") for fn in xrefs]
        try:
            logger.info("collecting crossreferences from %s files", len(xrefs))
            merge_sorted_h5_table(
                h5hs,
                table_path="/XRef",
                sort_columns=["EntryNr", "XRefSource", "XRefId", "Verification"],
                dupl_subset_columns=["EntryNr", "XRefSource", "XRefId"],
                storer_callback=storer.add_xrefs,
                enr_batch_size=100,
                chunksize=500_000,
            )
            logger.info("collecting EC annotations from %s files", len(xrefs))
            merge_sorted_h5_table(
                h5hs,
                table_path="/Annotations/EC",
                sort_columns=["EntryNr", "ECacc"],
                dupl_subset_columns=["EntryNr", "ECacc"],
                storer_callback=storer.add_ecs,
                enr_batch_size=1000,
                chunksize=500_000,
            )
        except Exception as e:
            logger.exception(f"Error while combining xrefs: {e}")
            raise
        finally:
            for h5h in h5hs:
                h5h.close()
