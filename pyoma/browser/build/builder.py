from __future__ import annotations

import collections
import csv
import errno
import fileinput
import gzip
import hashlib
import io
import itertools
import json
import multiprocessing as mp
import concurrent.futures
import operator
import os
import re
import math
import tempfile
import time
import codecs
from typing import Union, Optional, List, Iterable, Tuple, Literal

import numpy
import numpy.lib.recfunctions
import pandas
import tables
from PySAIS import sais

from tqdm import tqdm

from .. import suffixsearch
from .. import tablefmt
from ..KmerEncoder import KmerEncoder, DIGITS_AA
from ..OrthoXMLSplitter import OrthoXMLSplitter
from ..geneontology import GeneOntology, OntologyParser, FreqAwareGeneOntology
from ..homoeologs import HomeologsConfidenceCalculator
from ..synteny import SyntenyScorer
from .. import hoghelper
from .suffixarray_helper import build_filtered_sa, kmer_codes_for_positions
from ... import common, version
from ..convert import (
    uniq,
    hog_re,
    silentremove,
    gz_is_empty,
    load_tsv_to_numpy,
    read_vps_from_tsv,
    load_hogs_at_level,
    _load_taxonomy_without_ref_to_itselfs,
    compute_ortholog_types,
    get_or_create_tables_node,
    create_index_for_columns,
    create_fast_famhoglevel_lookup,
    create_and_store_fast_famhoglevel_lookup,
)
from ..exceptions import DBConsistencyError
from ..db import Taxonomy
from ..convert import DarwinExporter


# defining a handler the react on encoding errors while loading a file
def log_and_replace_error_handler(exception):
    # exception is a UnicodeDecodeError
    common.package_logger.error(
        f"Unicode decode error at byte {exception.start}: {exception.reason}. Replaceing it with '�'"
    )
    return ("�", exception.end)


codecs.register_error("logreplace", log_and_replace_error_handler)


class DataImportError(Exception):
    pass


class OmaGroupsProvider:
    def __init__(self, source):
        if source is not None:
            with open(source, "r") as f:
                self.data = json.load(f)
        else:
            self.data = None

    def get_oma_group(self, genome, nr):
        try:
            per_genome = self.data[str(genome)]
        except KeyError:
            return 0
        try:
            return per_genome[str(nr)]
        except KeyError:
            try:
                return per_genome[nr]
            except KeyError:
                return 0


class BufferedTableWriter:
    """
    Efficiently buffers and writes data to a PyTables table.
    Ensures buffer is a structured NumPy array for dtype safety and sorting.
    """

    def __init__(
        self,
        table,
        sort_key: Optional[List[str]] = None,
        buffer_size: int = 500_000,
        dtype: Optional[numpy.dtype] = None,
    ):
        self.table = table
        self.sort_key = sort_key
        self.buffer_size = buffer_size
        self.dtype = dtype or table.description._v_dtype
        self.buffer: List[numpy.array] = []
        self.buffer_cnt = 0
        self.total_written = 0

    def add(self, items: Union[Iterable[Tuple], numpy.ndarray, pandas.DataFrame]):
        """Add multiple items."""
        arr = self._to_recarray(items)
        if arr.size == 0:
            return

        self.buffer.append(arr)
        self.buffer_cnt += arr.size
        if self.buffer_cnt >= self.buffer_size:
            self.flush()

    def add_one(self, item: Tuple):
        """Add single row."""
        self.add([item])

    def _to_recarray(self, data):
        """Ensure data is structured NumPy array with proper dtype."""
        if isinstance(data, numpy.ndarray):
            return data.astype(self.dtype, copy=False)
        if isinstance(data, pandas.DataFrame):
            return data.to_records(index=False).astype(self.dtype, copy=False)
        # assume iterable of tuples
        return numpy.asarray(list(data), dtype=self.dtype)

    def flush(self):
        """Sort and write to table, then clear buffer."""
        if self.buffer_cnt == 0:
            return
        if len(self.buffer) == 1:
            buf = self.buffer[0]
        else:
            buf = numpy.concatenate(self.buffer)
        if self.sort_key:
            buf.sort(order=self.sort_key)
        self.table.append(buf)
        self.total_written += len(buf)
        self.buffer.clear()
        self.buffer_cnt = 0


class XrefFileWriter:
    """Handles one HDF5 file lifecycle and table management."""

    def __init__(self, filename, mode, index_cols=None, suffix_col=None, buffer_size=500_000):
        self.filename = filename
        self.mode = mode
        self.index_cols = index_cols or []
        self.suffix_col = suffix_col
        self.buffer_size = buffer_size
        self.h5 = None
        self.xref_writer = None
        self.ec_writer = None
        self.source_enum = None
        self.verify_enum = None

    def open(self):
        self.h5 = tables.open_file(
            self.filename, mode=self.mode, filters=tables.Filters(complevel=7, complib="blosc2", fletcher32=True)
        )
        if self.mode == "w":
            xref = self.h5.create_table("/", "XRef", tablefmt.XRefTable, expectedrows=10_000_000)
            ec = self.h5.create_table(
                "/Annotations", "EC", tablefmt.ECTable, expectedrows=1_000_000, createparents=True
            )
        else:
            xref = self.h5.root.XRef
            ec = self.h5.root.Annotations.EC

        self.xref_writer = BufferedTableWriter(
            xref, sort_key=["EntryNr", "XRefSource", "XRefId", "Verification"], buffer_size=self.buffer_size
        )
        self.ec_writer = BufferedTableWriter(ec, sort_key=["EntryNr", "ECacc"], buffer_size=self.buffer_size)

        self.source_enum = xref.get_enum("XRefSource")
        self.verify_enum = xref.get_enum("Verification")

    def write_xref(self, item):
        self.xref_writer.add_one(item)

    def write_xrefs(self, it):
        self.xref_writer.add(it)

    def write_ec(self, item):
        self.ec_writer.add_one(item)

    def write_ecs(self, it):
        self.ec_writer.add(it)

    def flush(self):
        self.xref_writer.flush()
        self.ec_writer.flush()

    @property
    def total_rows(self):
        return self.xref_writer.total_written + self.ec_writer.total_written

    def create_indexes(self):
        if self.index_cols:
            create_index_for_columns(self.h5.root.XRef, *self.index_cols)
            if "EntryNr" in self.index_cols:
                create_index_for_columns(self.h5.root.Annotations.EC, "EntryNr")
            if "XRefId" in self.index_cols:
                create_index_for_columns(self.h5.root.Annotations.EC, "ECacc")
        if self.suffix_col:
            suffixsearch.create_suffix_index(self.h5.root.XRef, self.suffix_col)

    def close(self):
        self.flush()
        if self.h5:
            self.h5.flush()
            self.create_indexes()
            self.h5.close()
            self.h5 = None


class XrefStorer:
    """
    Manages buffered writing of XRefs + ECs across multiple HDF5 files.
    Automatically rotates files after reaching max_rows_per_file.
    """

    def __init__(
        self,
        path: str,
        mode: Literal["a", "w", "r"] = "w",
        index_cols: Optional[List[str]] = None,
        suffix_col: Optional[str] = None,
        buffer_size: int = 500_000,
        multi_files: bool = False,
        max_rows_per_file: int = 10_000_000,
    ):
        self.path = path
        self.mode = mode
        if index_cols is not None:
            unknown_cols = set(index_cols) - set(tablefmt.XRefTable.columns.keys())
            if unknown_cols:
                raise ValueError("Unknown columns for indexing: {}".format(unknown_cols))
        self.index_cols = index_cols
        if suffix_col is not None:
            if suffix_col not in tablefmt.XRefTable.columns.keys():
                raise ValueError("Unknown columns building suffix index: {}".format(suffix_col))
        self.suffix_col = suffix_col
        self.buffer_size = buffer_size
        self.multi_files = multi_files
        self.max_rows_per_file = max_rows_per_file

        self._file_counter = 0 if multi_files else -1
        self._current_writer: Optional[XrefFileWriter] = None

    def __enter__(self):
        self._open_new_file_if_needed()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._current_writer:
            self._current_writer.close()

    def _build_filename(self):
        if self._file_counter < 0:
            return self.path
        base, ext = os.path.splitext(self.path)
        fname = f"{base}_{self._file_counter:03d}{ext}"
        self._file_counter += 1
        return fname

    def _open_new_file_if_needed(self):
        if not self._current_writer:
            fname = self._build_filename()
            self._current_writer = XrefFileWriter(
                fname, self.mode, self.index_cols, self.suffix_col, buffer_size=self.buffer_size
            )
            self._current_writer.open()
            return

        if self.multi_files and self._current_writer.total_rows >= self.max_rows_per_file:
            self._current_writer.flush()  # ensure all written
            self._current_writer.close()
            fname = self._build_filename()
            self._current_writer = XrefFileWriter(
                fname, self.mode, self.index_cols, self.suffix_col, buffer_size=self.buffer_size
            )
            self._current_writer.open()

    # ─────────────────────────────
    # Public interface
    # ─────────────────────────────
    def add_xref(self, enr: int, src: int, xref: str, verif: int, ident: float = 0.0):
        self._open_new_file_if_needed()
        item = (enr, src, xref.encode("utf-8"), verif, ident)
        self._current_writer.write_xref(item)

    def add_xrefs(self, it: Iterable[Tuple]):
        self._open_new_file_if_needed()
        self._current_writer.write_xrefs(it)

    def add_ec(self, enr: int, ec: str):
        self._open_new_file_if_needed()
        item = (enr, ec.encode("utf-8"))
        self._current_writer.write_ec(item)

    def add_ecs(self, it: Iterable[Tuple[int, bytes]]):
        self._open_new_file_if_needed()
        self._current_writer.write_ecs(it)

    def flush(self):
        if self._current_writer:
            self._current_writer.flush()

    def add_source_xref(self, enr: int, xref: str, typ: str):
        self._open_new_file_if_needed()
        src = (
            self._current_writer.source_enum["SourceID"]
            if typ == "id"
            else self._current_writer.source_enum["SourceAC"]
        )
        verif = self._current_writer.verify_enum["exact"]
        self.add_xref(enr, src, xref, verif, 1.0)


def load_homoeologs_from_tsv(genome, basedir: Optional[Union[str, os.PathLike]] = None):
    """load homoeologs from the tsv file"""
    if basedir is None or not os.path.isdir(basedir):
        common.package_logger.warning("No base directory for homoeologs passed. Won't load any homoeologs.")
        return numpy.empty(0, dtype=tables.dtype_from_descr(tablefmt.PairwiseRelationTable))
    fn = os.path.join(basedir, f"{genome['UniProtSpeciesCode'].decode()}.tsv.gz")
    off = genome["EntryOff"]
    if not os.path.exists(fn):
        common.package_logger.error("expected homoeologs file %s not found", fn)
        return numpy.empty(0, dtype=tables.dtype_from_descr(tablefmt.PairwiseRelationTable))
    return load_tsv_to_numpy((fn, off, off, False))


def identify_close_paralogs(df: pandas.DataFrame, join_threshold_mb=500) -> pandas.DataFrame:
    """identify the close paralogs (shared orthologs) pairs from a VPairs pandas.dataframe"""

    avg_group_size = df.groupby("EntryNr2").size().mean()
    num_groups = df["EntryNr2"].nunique()
    est_pairs = num_groups * (avg_group_size * (avg_group_size - 1)) / 2
    est_memory_bytes = est_pairs * 3 * 8
    est_memory_mb = est_memory_bytes / (1024**2)

    if est_memory_mb <= join_threshold_mb:
        # Join-based method
        common.package_logger.info(
            f"Using join-based method to identify close paralogs (estimated memory: {est_memory_mb:.1f} MB)"
        )
        df_with_index = df.set_index("EntryNr2")
        cp = df_with_index.join(df_with_index, rsuffix="_2")[["EntryNr1", "EntryNr1_2"]]
        cp = cp[cp["EntryNr1"] < cp["EntryNr1_2"]].drop_duplicates(ignore_index=True)
        cp = cp.rename(columns={"EntryNr1_2": "EntryNr2"})
    else:
        # Groupby-based method
        common.package_logger.info(
            f"⚠️ Using groupby-based method to identify close paralogs(estimated memory: {est_memory_mb:.1f} MB)"
        )

        # Group by EntryNr2
        grouped = df.groupby("EntryNr2")["EntryNr1"]

        # generate all pairs
        pairs = set()
        for entrynr2, entrynr1_group in grouped:
            entries = entrynr1_group.unique()
            if len(entries) > 1:
                pairs.update(pair for pair in itertools.combinations(sorted(entries), 2))
        cp = pandas.DataFrame(pairs, columns=["EntryNr1", "EntryNr2"])
    return cp.sort_values(by=["EntryNr1", "EntryNr2"], ignore_index=True)


class DBBuilder(DarwinExporter):
    def __init__(self, path, logger=None, mode=None, complib="zlib"):
        self.logger = logger if logger is not None else common.package_logger
        self._path = path
        self._complib = complib
        if mode is None:
            mode = "append" if os.path.exists(path) else "write"
        self._mode = mode

    def __enter__(self):
        compr = tables.Filters(complevel=7, complib=self._complib, fletcher32=False)
        self.h5 = tables.open_file(self._path, mode=self._mode[0], filters=compr)
        self.logger.info(f"opened {self._path} in {self._mode} mode, options {compr} ; pyoma {version()}")
        if self._mode == "write":
            self.h5.set_node_attr("/", "convertion_start", time.strftime("%c"))
            self.h5.set_node_attr("/", "pyoma_version", version())
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.h5.close()

    def call_darwin_export(self, func):
        raise NotImplementedError("Darwin export must not be called anymore")

    def get_version(
        self,
    ):
        """return version of the dataset.

        Default implementation searches for 'mname' in Matrix or matrix_stats.drw files.
        """
        return "Test"

    def add_taxonomy(self, tax_tsv):
        col_names = list(tablefmt.TaxonomyTable.columns)[:4]
        tax_data = pandas.read_csv(tax_tsv, sep="\t", names=col_names)
        dflt_cols = set(tablefmt.TaxonomyTable.columns) - set(tax_data.columns)
        for col in dflt_cols:
            tax_data[col] = tablefmt.TaxonomyTable.columns[col].dflt

        dt = {k: v.dtype for k, v in tablefmt.TaxonomyTable.columns.items()}
        taxtab = self.h5.create_table(
            "/",
            "Taxonomy",
            tablefmt.TaxonomyTable,
            obj=tax_data.to_records(index=False, column_dtypes=dt),
            expectedrows=len(tax_data),
        )
        create_index_for_columns(taxtab, "NCBITaxonId")

    def add_species_data(self, gs_tsv, taxid_updates=None):
        """parses a genome summary from the tsv file and adds it to the database"""

        def parse_as_date_column(val):
            if val == "":
                return 0
            for fmt in (
                "%b %d, %Y",
                "%B %d, %Y",
                "%d.%m.%Y",
                "%Y%m%d",
                "%Y-%m-%d",
                "%d-%m-%Y",
                "%d-%b-%Y",
            ):
                try:
                    date = time.strptime(val, fmt)
                    return time.mktime(date)
                except ValueError:
                    pass
            raise ValueError("Cannot parse date of '{}'".format(val))

        tax = Taxonomy(self.h5.get_node("/Taxonomy").read())
        taxid_order = {int(node["NCBITaxonId"]): i for i, (node, _) in enumerate(tax.traverse(strategy="postorder"))}

        tax_2_taxtable_row = {int(x["NCBITaxonId"]): x for x in tax.tax_table}

        def select_taxid_and_set_is_genome(row):
            taxid = row["NCBITaxonId"]
            if taxid not in tax_2_taxtable_row or not tax_2_taxtable_row[taxid]["IsGenome"]:
                taxid = row["GenomeId"]
            taxtabrow = tax_2_taxtable_row[taxid]
            return {"NCBITaxonId": taxid, "IsGenome": taxtabrow["IsGenome"], "SciName": taxtabrow["Name"]}

        data = pandas.read_csv(gs_tsv, sep="\t", dtype={"SciName": str})
        if taxid_updates is not None:
            data["NCBITaxonId"] = data["NCBITaxonId"].replace(taxid_updates)
        data[["NCBITaxonId", "taxid_is_genome", "SciName"]] = data.apply(
            select_taxid_and_set_is_genome, result_type="expand", axis=1
        )
        data["order"] = data["NCBITaxonId"].map(taxid_order)
        data.sort_values(by=["order", "GenomeId"], inplace=True)

        data.reset_index(drop=True, inplace=True)
        name2code = {str(row.Name): str(row.UniProtSpeciesCode) for row in data.itertuples(index=False)}

        cols = list(tablefmt.GenomeTable.columns)
        dflt_cols = set(cols) - set(data.columns)
        for col in dflt_cols:
            data[col] = tablefmt.GenomeTable.columns[col].dflt

        # Build EntryOff after sorting genomes
        data.loc[0, "EntryOff"] = 0
        for i in range(len(data) - 1):
            data.loc[i + 1, "EntryOff"] = data.loc[i, "EntryOff"] + data.loc[i, "TotEntries"]

        gs = data[cols]
        for col, typeinfo in tablefmt.GenomeTable.columns.items():
            if typeinfo.kind == "time":
                gs.loc[:, col] = gs.loc[:, col].apply(parse_as_date_column)
            elif typeinfo.kind == "string":
                gs.loc[:, col] = gs.loc[:, col].fillna("")
        dt = {k: v.dtype for k, v in tablefmt.GenomeTable.columns.items()}
        gstab = self.h5.create_table(
            "/", "Genome", tablefmt.GenomeTable, obj=gs.to_records(index=False, column_dtypes=dt), expectedrows=len(gs)
        )
        create_index_for_columns(gstab, "NCBITaxonId", "UniProtSpeciesCode", "EntryOff")
        return name2code

    def add_orthologs(
        self,
        basedir: Optional[Union[str, os.PathLike]],
        genomes: tables.Table,
        homoeologs_base: Optional[Union[str, os.PathLike]],
    ):
        genome_offs = genomes.col("EntryOff")
        anygenome = genomes[0]["UniProtSpeciesCode"].decode()
        if basedir is None:
            self.logger.warning(
                "No base directory for pairwise orthologs passed. Will initialize with empty VPair tables"
            )
        else:
            testdir = os.path.join(basedir, anygenome)
            if not (os.path.isdir(testdir) and any(map(lambda x: x.endswith(".orth.txt.gz"), os.listdir(testdir)))):
                raise RuntimeError(f"{basedir} does not contain ortholog files")
            self.logger.info("using %s as base dir for pairwise orthology", basedir)

        for gs in genomes.iterrows():
            genome = gs["UniProtSpeciesCode"].decode()
            rel_node_for_genome = self._get_or_create_node(f"/PairwiseRelation/{genome}")
            if "VPairs" not in rel_node_for_genome:
                if basedir is None:
                    data = []
                else:
                    data = read_vps_from_tsv(
                        genomes, genome.encode("utf-8"), basedir=basedir, check_exist_and_swap=True
                    )
                vp_tab = self.h5.create_table(
                    rel_node_for_genome,
                    "VPairs",
                    tablefmt.PairwiseRelationTable,
                    expectedrows=len(data),
                )
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, vp_tab)
                if numpy.any(data["RelType"] >= tablefmt.PairwiseRelationTable.columns.get("RelType").enum["n/a"]):
                    compute_ortholog_types(data, genome_offs)
                self._write_to_table(vp_tab, data)
                create_index_for_columns(vp_tab, "EntryNr1")
            if "within" not in rel_node_for_genome:
                df = pandas.DataFrame(self.h5.get_node(rel_node_for_genome, "VPairs").read())
                df_with_ss_paralogs = df.loc[df["RelType"] > 1, ["EntryNr1", "EntryNr2"]]
                cp = identify_close_paralogs(df_with_ss_paralogs)
                cp["RelType"] = tablefmt.PairwiseRelationTable.columns.get("RelType").enum["close paralog"]
                cols = list(tablefmt.PairwiseRelationTable.columns)
                dflt_cols = set(cols) - set(cp.columns)
                for col in dflt_cols:
                    cp[col] = tablefmt.PairwiseRelationTable.columns[col].dflt
                cp = cp[cols]
                if gs["IsPolyploid"] and homoeologs_base is not None:
                    hp = pandas.DataFrame(load_homoeologs_from_tsv(gs, basedir=homoeologs_base))
                    dflt_cols = set(cols) - set(hp.columns)
                    for col in dflt_cols:
                        hp[col] = tablefmt.PairwiseRelationTable.columns[col].dflt
                    hp = hp.set_index(["EntryNr1", "EntryNr2"])
                    cp = cp.set_index(["EntryNr1", "EntryNr2"])
                    cp.update(hp)
                    cp = cp.reset_index(drop=False)
                    cp = cp[cols]

                dt = {k: v.dtype for k, v in tablefmt.PairwiseRelationTable.columns.items()}
                within_tab = self.h5.create_table(
                    rel_node_for_genome,
                    "within",
                    tablefmt.PairwiseRelationTable,
                    obj=cp.to_records(index=False, column_dtypes=dt),
                    expectedrows=len(cp),
                )
                create_index_for_columns(within_tab, "EntryNr1")

    def _add_desc(self, desc, row, array):
        row["DescriptionOffset"] = len(array)
        row["DescriptionLength"] = len(desc)
        desc_arr = numpy.ndarray((len(desc),), buffer=desc.encode("utf-8"), dtype=tables.StringAtom(1))
        array.append(desc_arr)

    def add_proteins(self, code_to_file, oma_group_provider, xref_collector):
        gs_node = self.h5.get_node("/Genome")
        if len(code_to_file) < len(gs_node):
            raise ValueError(
                f"nr of json files does not match number of genomes: " f"{len(code_to_file)} vs {len(gs_node)}"
            )
        nr_prot = sum(gs_node.cols.TotEntries)
        nr_aa = sum(gs_node.cols.TotAA)
        prot_grp = self._get_or_create_node("/Protein", "Root node for protein (oma entries) information")
        prot_tab = self.h5.create_table(prot_grp, "Entries", tablefmt.ProteinTable, expectedrows=nr_prot)
        seq_arr = self.h5.create_earray(
            prot_grp,
            "SequenceBuffer",
            tables.StringAtom(1),
            (0,),
            "concatenated protein sequences",
            expectedrows=nr_aa + nr_prot,
        )
        cdna_arr = self.h5.create_earray(
            prot_grp,
            "CDNABuffer",
            tables.StringAtom(1),
            (0,),
            "concatenated cDNA sequences",
            expectedrows=3 * nr_aa + nr_prot,
        )
        desc_arr = self.h5.create_earray(
            prot_grp,
            "DescriptionBuffer",
            tables.StringAtom(1),
            (0,),
            "concatenated protein descriptions",
            expectedrows=100 * nr_prot,
        )

        for gs in gs_node.iterrows():
            genome = gs["UniProtSpeciesCode"].decode()
            with open(code_to_file[genome], encoding="utf-8", errors="logreplace") as fd:
                data = json.load(fd)
            if len(data["seqs"]) != gs["TotEntries"]:
                raise DataImportError(
                    f"number of entries ({len(data['seqs'])}) does not match number "
                    f"of seqs ({gs['TotEntries']}) for {genome}"
                )

            loc_tab = self.h5.create_table(
                "/Protein/Locus",
                genome,
                tablefmt.LocusTable,
                createparents=True,
                expectedrows=sum(len(z) for z in data["locs"]),
            )

            cnt_missmatch_locus = 0
            for nr in range(gs["TotEntries"]):
                e_nr = gs["EntryOff"] + nr + 1
                prot_tab.row["EntryNr"] = e_nr
                prot_tab.row["OmaGroup"] = oma_group_provider.get_oma_group(genome, nr + 1)

                self._add_sequence(data["seqs"][nr], prot_tab.row, seq_arr)
                self._add_sequence(data["cdna"][nr], prot_tab.row, cdna_arr, "CDNA")
                self._add_desc(data["de"][nr], prot_tab.row, desc_arr)

                prot_tab.row["Chromosome"] = data["chrs"][nr]
                prot_tab.row["OmaHOG"] = b""  # will be assigned later
                prot_tab.row["CanonicalId"] = data["acs"][nr][0].encode("utf-8")
                if xref_collector is not None:
                    for xref in data["acs"][nr]:
                        xref_collector.add_source_xref(e_nr, xref, "ac")
                    xref_collector.add_source_xref(e_nr, data["ids"][nr], "id")

                # prot_tab.row["AltSpliceVariant"] = data["alts"][nr]
                # if prot_tab.row["AltSpliceVariant"] == 0 or prot_tab.row["AltSpliceVariant"] == prot_tab.row["EntryNr"]:
                #     cnt_genes += 1  # main isoforms of gene

                locus_tab = numpy.array([(e_nr, *row) for row in data["locs"][nr]], dtype=loc_tab.dtype)
                loc_tab.append(locus_tab)
                len_cds = sum(z["End"] - z["Start"] + 1 for z in locus_tab)
                if len_cds != prot_tab.row["CDNABufferLength"] - 1:
                    if cnt_missmatch_locus < 10:
                        self.logger.debug(
                            f"Sum of exon lengths differs cDNA sequence {genome}{nr+1:05d} ({len(locus_tab)} exons): "
                            f"{len_cds} vs {prot_tab.row['CDNABufferLength'] - 1}"
                        )
                    cnt_missmatch_locus += 1
                if len(locus_tab) > 0:
                    prot_tab.row["LocusStart"] = locus_tab["Start"].min()
                    prot_tab.row["LocusEnd"] = locus_tab["End"].max()
                    prot_tab.row["LocusStrand"] = locus_tab[0]["Strand"]
                else:
                    prot_tab.row["LocusStart"] = 0
                    prot_tab.row["LocusEnd"] = 0
                    prot_tab.row["LocusStrand"] = 1
                if gs["IsPolyploid"]:
                    prot_tab.row["SubGenome"] = data["subgenome"][nr].encode("ascii")
                prot_tab.row.append()
            prot_tab.flush()
            seq_arr.flush()
            # gs["TotGenes"] = cnt_genes
            # gs.update()
            if cnt_missmatch_locus > 0:
                self.logger.warning(
                    "[%s]: %d miss-matches in exon-lengths compared to locus info", genome, cnt_missmatch_locus
                )
            for n in (prot_tab, seq_arr, loc_tab):
                if n.size_in_memory != 0:
                    self.logger.info(
                        "worte %s: compression ratio %3f%%" % (n._v_pathname, 100 * n.size_on_disk / n.size_in_memory)
                    )
        create_index_for_columns(prot_tab, "EntryNr", "MD5ProteinHash", "OmaGroup")

    def add_sequence_index(self, seqs: bytes, nr_entries: int, k: int = 6):
        """compute the suffix array and Kmer lookup index and store in the hdf5 under /Protein
        :param seqs: concatenated sequences, delimitted between entries with a space
        :param nr_entries: number of entries in the database
        :param k: kmer size used for KmerIndex.
        """
        # Compute & save the suffix array to DB.
        self.logger.info("computing suffix array")
        sa = sais(seqs)
        sa[:nr_entries].sort()  # Sort delimiters by position.

        self.logger.info(f"storing suffix array into database ({len(sa)} positions)")
        sa_h5 = self.h5.create_carray(
            "/Protein",
            createparents=True,
            name="SequenceIndex",
            title="concatenated protein sequences suffix array",
            obj=sa,
        )
        delimiters_sorted = sa[:nr_entries]

        # store dtype and atom type to fit size of actual entry number and suffix array
        dtype_sa = sa.dtype
        dtype_enr = numpy.uint32 if (nr_entries < numpy.iinfo(numpy.uint32).max) else numpy.uint64

        sa = None  # free memory

        n_seq = len(seqs)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = os.path.join(tmpdir, "sa_helper.h5")
            h5_tmp = tables.open_file(
                tmp_path, mode="w", filters=tables.Filters(complevel=3, complib="blosc2", bitshuffle=True)
            )

            # -----------------------------
            # Phase 2: Stream-filter SA → SA_filtered, SA_origpos, and build IndexInSAOrder
            # -----------------------------
            sa_f, sa_origpos, idx_saorder = build_filtered_sa(
                sa_h5, delimiters_sorted, n_seq, nr_entries, k=k, dtype_sa=dtype_sa, dtype_enr=dtype_enr, h5_out=h5_tmp
            )

            # -----------------------------
            # Phase 3: Build k-mer lookup by scanning SA_filtered
            # -----------------------------
            kmers = KmerEncoder(k, is_protein=True)
            kmer_lookup_arr = self.h5.create_vlarray(
                "/Protein",
                name="KmerLookup",
                atom=tables.Atom.from_dtype(numpy.dtype(dtype_enr)),
                title="kmer entry lookup table",
                expectedrows=len(kmers),
            )
            self.h5.set_node_attr(kmer_lookup_arr, "k", k)

            # vectorized preparation
            seqs_np = numpy.frombuffer(seqs, dtype=numpy.uint8)
            map256 = numpy.full(256, 255, dtype=numpy.uint8)
            for i, aa in enumerate(DIGITS_AA):  # or DIGITS_DNA if not protein
                map256[ord(aa)] = i
            alphabet_size = len(DIGITS_AA)

            # helper to grow VLArray up to a code index
            def _ensure_rows(vl, upto: int):
                need = upto - len(vl)
                for _ in range(need):
                    vl.append([])

            def _emit_run(code_int: int, arr: numpy.ndarray):
                _ensure_rows(kmer_lookup_arr, code_int)
                kmer_lookup_arr.append(arr)

            chunksize = sa_h5.chunkshape[0] * (
                1
                if len(sa_h5) // sa_h5.chunkshape[0] < 64_000
                else math.ceil(len(sa_h5) / sa_h5.chunkshape[0] / 64_000)
            )

            # Now find the split points and construct lookup ragged array.
            L = int(len(sa_f)) - k
            t = tqdm(total=max(L, 0), desc="Building Kmer lookup")
            ii, tot_kmers = 0, len(kmers)
            chunk_keys: list[int] = []  # chunk keys for sa-lookup-table
            chunk_pos: list[int] = []  # positions of chunk in sa-lookup-table

            # carry-over for runs and for a deferred boundary that fell inside a run crossing batches
            prev_code: int | None = None
            prev_entries: list[numpy.ndarray] = []  # across batch kmers
            pending_cut: bool = False

            if L > 0:
                first_code = int(kmer_codes_for_positions(sa_f[0:1], k, seqs_np, dtype_sa, map256, alphabet_size)[0])
                chunk_pos.append(int(sa_origpos[0]))
                chunk_keys.append(first_code)
            next_target = chunk_pos[-1] + int(chunksize) if chunk_pos else int(chunksize)

            batch_size = 2**16
            while ii < L:
                jj = min(ii + batch_size, L)
                P = sa_f[ii:jj]  # start positions
                codes, good = kmer_codes_for_positions(P, k, seqs_np, dtype_sa, map256, alphabet_size)
                entries = idx_saorder[ii:jj]
                origpos = sa_origpos[ii:jj]
                # run boundaries inside the batch

                if jj - ii > 1:
                    change = numpy.nonzero(codes[1:] != codes[:-1])[0] + 1
                    run_starts = numpy.concatenate(([0], change))
                    run_ends = numpy.concatenate((change, [codes.size]))
                else:
                    run_starts = numpy.array([0], dtype=numpy.int64)
                    run_ends = numpy.array([1], dtype=numpy.int64)

                # is the first local run a continuation of previous batch?
                first_is_cont = prev_code is not None and int(codes[0]) == prev_code

                for a, b in zip(run_starts, run_ends):
                    code_int = int(codes[a])

                    # entries slice for this run (for VLArray emission)
                    entries_slice = entries[a:b].astype(dtype_enr, copy=False)

                    # ---- boundary scheduler (align to run starts) ----
                    # We may need to place 0..N boundaries while scanning forward.
                    # A boundary belongs at the *start of a run*.
                    start_pos = int(origpos[a])
                    end_pos = int(origpos[b - 1])

                    # 1) If we were waiting for the *next* run start (because the target fell inside
                    #    a run that continued past the previous batch), place it now at this run start.
                    if pending_cut:
                        chunk_pos.append(start_pos)
                        chunk_keys.append(code_int)
                        next_target = chunk_pos[-1] + int(chunksize)
                        pending_cut = False  # satisfied by this run start

                    # 2) Place as many boundaries as targets we cross while scanning.
                    #    A boundary can be placed at:
                    #      - this run's start (if target <= start_pos and this is a *true* run start),
                    #      - otherwise at the *next* run start (i.e., end of this run).
                    while next_target <= end_pos:
                        this_is_true_start = not (a == 0 and first_is_cont)
                        if this_is_true_start and next_target <= start_pos:
                            # cut exactly at this run start
                            chunk_pos.append(start_pos)
                            chunk_keys.append(code_int)
                            next_target = chunk_pos[-1] + int(chunksize)
                            # only one boundary per run start; further cuts must wait for later runs
                            break
                        else:
                            # target lies *inside* this run (or at a fake start due to continuation)
                            # → align cut to the *next* run start (code change).
                            if b < len(codes):
                                next_start_pos = int(origpos[b])  # start of next run in this batch
                                next_code_int = int(codes[b])
                                chunk_pos.append(next_start_pos)
                                chunk_keys.append(next_code_int)
                                next_target = chunk_pos[-1] + int(chunksize)
                                # we may still have more targets to satisfy inside subsequent runs,
                                # but not inside the current run (we aligned to its end), so break.
                                break
                            else:
                                # the next run start is in the *next batch* → defer
                                pending_cut = True
                                # do not append now; it will be appended at the next run start we see
                                break

                    # ---- emit VLArray data, merging with a possible carry-over run ----
                    if prev_code is None:
                        prev_code = code_int
                        prev_entries = [entries_slice]
                    elif code_int == prev_code:
                        prev_entries.append(entries_slice)
                    else:
                        _emit_run(prev_code, numpy.concatenate(prev_entries))
                        prev_code = code_int
                        prev_entries = [entries_slice]

                t.update(jj - ii)
                ii = jj

            # flush final run
            if prev_code is not None and prev_entries:
                _emit_run(prev_code, numpy.concatenate(prev_entries))

            # pad remaining empty rows
            for _ in range(tot_kmers - len(kmer_lookup_arr)):
                kmer_lookup_arr.append([])
            kmer_lookup_arr.flush()
            t.close()
            h5_tmp.close()

        self.logger.info("storing suffix array lookup index into database")
        chunk_pos.append(len(seqs))
        chunk_keys.append(tot_kmers)
        self.h5.create_carray("/Protein", name="SuffixArrayIndexKeys", obj=chunk_keys, title="suffix array keys")
        self.h5.create_carray("/Protein", name="SuffixArrayIndexPos", obj=chunk_pos, title="suffix array positions")

    def add_protein_hog_ids(self, hog_ids: numpy.array) -> None:
        entries_tab = self.h5.get_node("/Protein/Entries")
        assert len(hog_ids) == len(entries_tab)
        entries_tab.modify_column(0, len(entries_tab), 1, column=hog_ids, colname="OmaHOG")
        create_index_for_columns(entries_tab, "OmaHOG")

    def identify_and_store_splice_variants(self, splice_json):
        with open(splice_json, "rt") as f:
            splice_info = json.load(f)
        gs = {
            row["UniProtSpeciesCode"].decode(): slice(int(row["EntryOff"]), int(row["EntryOff"] + row["TotEntries"]), 1)
            for row in self.h5.get_node("/Genome")
        }
        entry_tab = self.h5.get_node("/Protein/Entries")
        alt_splice = numpy.zeros((len(entry_tab),), dtype=entry_tab.cols.AltSpliceVariant.dtype)
        for sp in splice_info:
            if sp not in gs:
                continue
            self._identify_main_variants(
                splice_groups=splice_info[sp],
                splice_arr=alt_splice,
                entries=entry_tab.read(start=gs[sp].start, stop=gs[sp].stop),
                offset=gs[sp].start,
                vp_tab=self.h5.get_node(f"/PairwiseRelation/{sp}/VPairs"),
            )
        entry_tab.modify_column(0, len(entry_tab), 1, column=alt_splice, colname="AltSpliceVariant")
        entry_tab.flush()
        self.update_nr_genes()

    def _identify_main_variants(self, splice_groups, splice_arr, entries, offset, vp_tab):
        for grp in splice_groups:
            idx = numpy.array(grp, dtype="i4") - 1
            ent = entries[idx]
            og = numpy.nonzero(ent["OmaGroup"])[0]
            if len(og) > 1:
                raise DBConsistencyError("Several splice variants in OMA Groups", ent)
            if len(og) == 1:
                splice_arr[idx + offset] = ent["EntryNr"][og[0]]
                continue

            hog = numpy.nonzero(ent["OmaHOG"])[0]
            if len(hog) > 1:
                raise DBConsistencyError("Several splice variants in HOGs", ent)
            if len(hog) == 1:
                splice_arr[idx + offset] = ent["EntryNr"][hog[0]]
                continue

            nr_vps = numpy.fromiter(
                map(lambda enr: common.count_elements(vp_tab.where("EntryNr1 == enr")), ent["EntryNr"]), dtype="i4"
            )
            vp = numpy.nonzero(nr_vps)[0]
            if len(vp) > 1:
                raise DBConsistencyError("Several splice variants contain pairwise orthologs", ent, nr_vps)
            if len(vp) == 1:
                splice_arr[idx + offset] = ent["EntryNr"][vp[0]]

            # no orthologs for any variant. choose the longest variant as main one.
            splice_arr[idx + offset] = ent["EntryNr"][numpy.argmax(ent["SeqBufferLength"])]

    def update_nr_genes(self):
        etab = self.h5.get_node("/Protein/Entries")
        for row in self.h5.get_node("/Genome"):
            rng = slice(row["EntryOff"], row["EntryOff"] + row["TotEntries"], 1)
            row["TotGenes"] = common.count_elements(
                etab.where("(EntryNr == AltSpliceVariant) | (AltSpliceVariant == 0)", start=rng.start, stop=rng.stop)
            )
            row.update()
        etab.flush()
