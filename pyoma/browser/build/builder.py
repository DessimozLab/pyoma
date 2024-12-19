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
import time
from typing import Union, Optional, List, Iterable, Tuple

import numpy
import numpy.lib.recfunctions
import pandas
import tables
from PySAIS import sais

from tqdm import tqdm

from .. import suffixsearch
from .. import tablefmt
from ..KmerEncoder import KmerEncoder
from ..OrthoXMLSplitter import OrthoXMLSplitter
from ..geneontology import GeneOntology, OntologyParser, FreqAwareGeneOntology
from ..homoeologs import HomeologsConfidenceCalculator
from ..synteny import SyntenyScorer
from .. import hoghelper
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
            return self.data[str(genome)][str(nr)]
        except KeyError:
            return 0


class XrefStorer:
    def __init__(self, path, mode: str = "w", index_cols: Optional[List] = None, suffix_col: Optional[str] = None):
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

    def __enter__(self):
        self.h5 = tables.open_file(
            self.path, mode=self.mode, filters=tables.Filters(complevel=5, complib="blosc2", fletcher32=True)
        )
        if self.mode == "w":
            self.xref = self.h5.create_table("/", "XRef", tablefmt.XRefTable, expectedrows=1e7)
            self.ec = self.h5.create_table("/Annotations", "EC", tablefmt.ECTable, expectedrows=1e6, createparents=True)
        self.source_enum = self.xref.get_enum("XRefSource")
        self.verify_enum = self.xref.get_enum("Verification")
        self._buffer = []
        self._ecbuffer = []
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.flush()
        self.h5.flush()
        if self.index_cols is not None:
            create_index_for_columns(self.xref, *self.index_cols)
            if "EntryNr" in self.index_cols:
                create_index_for_columns(self.ec, "EntryNr")
            if "XRefId" in self.index_cols:
                create_index_for_columns(self.ec, "ECacc")
        if self.suffix_col is not None:
            suffixsearch.create_suffix_index(self.xref, self.suffix_col)
        self.h5.close()

    def flush(self):
        if len(self._buffer) > 0:
            self.xref.append(self._buffer)
        if len(self._ecbuffer) > 0:
            self.ec.append(self._ecbuffer)
        self._buffer = []
        self._ecbuffer = []

    def add_xrefs(self, it: Iterable[Tuple]):
        """adds a bunch of xrefs from it. The tuple must be in the correct format (no checks performed)"""
        self._buffer.extend(it)
        if len(self._buffer) > 500_000:
            self.flush()

    def add_xref(self, enr: int, src: int, xref: str, verif: int, ident: float = 0):
        """adds an xref entry. src and verif need to be already mapped to their numeric enum value"""
        self._buffer.append((enr, src, xref.encode("utf-8"), verif, ident))
        if len(self._buffer) > 500_000:
            self.flush()

    def add_ecs(self, it: Iterable[Tuple[int, bytes]]):
        self._ecbuffer.extend(it)
        if len(self._ecbuffer) > 50_000:
            self.flush()

    def add_ec(self, enr, ec):
        self._ecbuffer.append((enr, ec.encode("utf-8")))
        if len(self._ecbuffer) > 50_000:
            self.flush()

    def add_source_xref(self, enr: int, xref: str, typ: str):
        """Adds a source xref entry. type needs to be either 'id' or 'ac'."""
        src = self.source_enum["SourceID"] if typ == "id" else self.source_enum["SourceAC"]
        self.add_xref(enr, src, xref, self.verify_enum["exact"], 1)


class DBBuilder(DarwinExporter):
    def __init__(self, path, logger=None, mode=None, complib="zlib"):
        self.logger = logger if logger is not None else common.package_logger
        self._path = path
        self._complib = complib
        if mode is None:
            mode = "append" if os.path.exists(path) else "write"
        self._mode = mode

    def __enter__(self):
        compr = tables.Filters(complevel=6, complib=self._complib, fletcher32=False)
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

    def add_species_data(self, gs_tsv):
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

        data = pandas.read_csv(gs_tsv, sep="\t")
        data["order"] = data["NCBITaxonId"].map(taxid_order)
        data.sort_values(by=["order", "GenomeId"], inplace=True)
        data["NCBITaxonId"] = numpy.where(
            (data["NCBITaxonId"].duplicated(keep=False)), data["GenomeId"], data["NCBITaxonId"]
        )
        data.reset_index(drop=True, inplace=True)

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

    def add_orthologs(self, basedir: Union[str, os.PathLike], genomes: tables.Table):
        genome_offs = genomes.col("EntryOff")
        anygenome = genomes[0]["UniProtSpeciesCode"].decode()
        testdir = os.path.join(basedir, anygenome)
        if not (os.path.isdir(testdir) and any(map(lambda x: x.endswith(".orth.txt.gz"), os.listdir(testdir)))):
            raise RuntimeError(f"{basedir} does not contain ortholog files")

        self.logger.info("using %s as base dir for pairwise orthology", basedir)
        for gs in genomes.iterrows():
            genome = gs["UniProtSpeciesCode"].decode()
            rel_node_for_genome = self._get_or_create_node(f"/PairwiseRelation/{genome}")
            if "VPairs" not in rel_node_for_genome:
                data = read_vps_from_tsv(genomes, genome.encode("utf-8"), basedir=basedir, check_exist_and_swap=True)
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
                df_with_ss_paralogs = df.loc[df["RelType"] > 1, ["EntryNr1", "EntryNr2"]].set_index("EntryNr2")
                cp = df_with_ss_paralogs.join(df_with_ss_paralogs, rsuffix="_2")
                cp = cp[cp["EntryNr1"] < cp["EntryNr1_2"]].reset_index()
                cp["EntryNr2"] = cp["EntryNr1_2"]
                cp["RelType"] = tablefmt.PairwiseRelationTable.columns.get("RelType").enum["close paralog"]
                cols = list(tablefmt.PairwiseRelationTable.columns)
                dflt_cols = set(cols) - set(cp.columns)
                for col in dflt_cols:
                    cp[col] = tablefmt.PairwiseRelationTable.columns[col].dflt
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

    def add_proteins(self, genome_files, oma_group_provider, xref_collector):
        code_to_file = {os.path.basename(f).split(".")[0]: f for f in genome_files}
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
            with open(code_to_file[genome], "r") as fd:
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
                prot_tab.row["LocusStart"] = locus_tab["Start"].min()
                prot_tab.row["LocusEnd"] = locus_tab["End"].max()
                prot_tab.row["LocusStrand"] = locus_tab[0]["Strand"]
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
        sa = sais(seqs)
        sa[:nr_entries].sort()  # Sort delimiters by position.
        self.h5.create_carray(
            "/Protein",
            createparents=True,
            name="SequenceIndex",
            title="concatenated protein sequences suffix array",
            obj=sa,
        )

        # Create lookup table for fa2go
        dtype = numpy.uint32 if (nr_entries < numpy.iinfo(numpy.uint32).max) else numpy.uint64
        idx = numpy.zeros(sa.shape, dtype=dtype)
        mask = numpy.zeros(sa.shape, dtype=bool)

        # Compute mask and entry index for sequence buff
        for i in range(nr_entries):
            s = (sa[i - 1] if i > 0 else -1) + 1
            e = sa[i] + 1
            idx[s:e] = i + 1
            mask[(e - k) : e] = True  # (k-1) invalid and delim.

        # Mask off those we don't want...
        sa = sa[~mask[sa]]

        # Reorder the necessary elements of entry index
        idx = idx[sa]

        # Initialise lookup array
        atom = tables.UInt32Atom if dtype is numpy.uint32 else tables.UInt64Atom
        kmers = KmerEncoder(k, is_protein=True)
        kmer_lookup_arr = self.h5.create_vlarray(
            "/Protein",
            name="KmerLookup",
            atom=atom(shape=()),
            title="kmer entry lookup table",
            expectedrows=len(kmers),
        )
        self.h5.set_node_attr(kmer_lookup_arr, "k", k)

        # Now find the split points and construct lookup ragged array.
        ii = 0
        for kk in tqdm(range(len(kmers)), desc="Constructing kmer lookup"):
            kmer = kmers.encode(kk)
            if (ii < len(sa)) and (seqs[sa[ii] : (sa[ii] + k)] == kmer):
                jj = ii + 1
                while (jj < len(sa)) and (seqs[sa[jj] : (sa[jj] + k)] == kmer):
                    jj += 1
                kmer_lookup_arr.append(idx[ii:jj])
                # New start
                ii = jj
            else:
                # End or not found
                kmer_lookup_arr.append([])
        kmer_lookup_arr.flush()

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
