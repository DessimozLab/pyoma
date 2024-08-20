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

import familyanalyzer
import lxml.html
import numpy
import numpy.lib.recfunctions
import pandas
import tables
from PySAIS import sais
from future.standard_library import hooks
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
from ..convert import DarwinExporter


class DataImportError(Exception):
    pass


class DBBuilder(DarwinExporter):
    def __init__(self, path, logger=None, mode=None):
        self.logger = logger if logger is not None else common.package_logger
        self._path = path
        if mode is None:
            mode = "append" if os.path.exists(path) else "write"
        self._mode = mode

    def __enter__(self):
        compr = tables.Filters(complevel=6, complib="zlib", fletcher32=True)
        self.h5 = tables.open_file(self._path, mode=self._mode[0], filters=compr)
        self.logger.info(f"opened {self._path} in {self._mode} mode, options {compr} ; pyoma {version()}")
        if self._mode == "write":
            self.h5.set_node_attr("/", "convertion_start", time.strftime("%c"))
            self.h5.set_node_attr("/", "pyoma_version", version())

    def call_darwin_export(self, func):
        raise NotImplementedError("Darwin export must not be called anymore")

    def get_version(
        self,
    ):
        """return version of the dataset.

        Default implementation searches for 'mname' in Matrix or matrix_stats.drw files.
        """
        return "Test"

    def add_species_data(self, gs_tsv, tax_tsv):
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

        data = pandas.read_csv(gs_tsv, sep="\t")
        cols = list(tablefmt.GenomeTable.columns)
        dflt_cols = set(cols) - set(data.columns)
        for col in dflt_cols:
            data[col] = tablefmt.GenomeTable.columns[col].dflt

        gs = data[cols]
        for col, typeinfo in tablefmt.GenomeTable.columns.items():
            if typeinfo.kind == "time":
                gs[col] = gs[col].apply(parse_as_date_column)
        dt = {k: v.dtype for k, v in tablefmt.GenomeTable.columns.items()}
        gstab = self.h5.create_table(
            "/", "Genome", tablefmt.GenomeTable, obj=gs.to_records(index=False, column_dtypes=dt), expectedrows=len(gs)
        )
        create_index_for_columns(gstab, "NCBITaxonId", "UniProtSpeciesCode", "EntryOff")

        col_names = list(tablefmt.TaxonomyTable.columns)[:3]
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

    def add_orthologs(self, basedir):
        genome_offs = self.h5.root.Genome.col("EntryOff")
        anygenome = self.h5.root.Genome[0]["UniProtSpeciesCode"].decode()
        testdir = os.path.join(basedir, anygenome)
        if not (os.path.isdir(testdir) and any(map(lambda x: x.endswith(".orth.txt.gz"), os.listdir(testdir)))):
            raise RuntimeError(f"{basedir} does not contain ortholog files")

        self.logger.info("using %s as base dir for pairwise orthology", basedir)
        for gs in self.h5.root.Genome.iterrows():
            genome = gs["UniProtSpeciesCode"].decode()
            rel_node_for_genome = self._get_or_create_node("/PairwiseRelation/{}".format(genome))
            if "VPairs" not in rel_node_for_genome:
                data = read_vps_from_tsv(self.h5.root.Genome, genome.encode("utf-8"), basedir=basedir)
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

    def _add_sequence(self, sequence, row, sequence_array, off, typ="Seq"):
        # add ' ' after each sequence (Ascii is smaller than
        # any AA, allows to build PAT array with split between
        # sequences.
        seqLen = len(sequence) + 1
        row[typ + "BufferOffset"] = off
        row[typ + "BufferLength"] = seqLen
        if typ == "CDNA":
            sequence = sequence.replace("X", "N")
        seqNumpyObj = numpy.ndarray(
            (seqLen,),
            buffer=(sequence + " ").encode("utf-8"),
            dtype=tables.StringAtom(1),
        )
        sequence_array.append(seqNumpyObj)
        if typ == "Seq":
            row["MD5ProteinHash"] = hashlib.md5(sequence.encode("utf-8")).hexdigest()
        return seqLen

    def add_proteins(self):
        gsNode = self.h5.get_node("/Genome")
        nrProt = sum(gsNode.cols.TotEntries)
        nrAA = sum(gsNode.cols.TotAA)
        protGrp = self._get_or_create_node("/Protein", "Root node for protein (oma entries) information")
        protTab = self.h5.create_table(protGrp, "Entries", tablefmt.ProteinTable, expectedrows=nrProt)
        seqArr = self.h5.create_earray(
            protGrp,
            "SequenceBuffer",
            tables.StringAtom(1),
            (0,),
            "concatenated protein sequences",
            expectedrows=nrAA + nrProt,
        )
        cdnaArr = self.h5.create_earray(
            protGrp,
            "CDNABuffer",
            tables.StringAtom(1),
            (0,),
            "concatenated cDNA sequences",
            expectedrows=3 * nrAA + nrProt,
        )
        seqOff = cdnaOff = 0
        loc_parser = locus_parser.LocusParser()
        for gs in gsNode.iterrows():
            genome = gs["UniProtSpeciesCode"].decode()
            cache_file = os.path.join(
                os.getenv("DARWIN_NETWORK_SCRATCH_PATH", ""),
                "pyoma",
                "prots",
                "{}.json".format(genome),
            )
            if os.path.exists(cache_file):
                with open(cache_file, "r") as fd:
                    data = json.load(fd)
            else:
                data = self.call_darwin_export("GetProteinsForGenome({})".format(genome))

            if len(data["seqs"]) != gs["TotEntries"]:
                raise DataImportError(
                    "number of entries ({:d}) does "
                    "not match number of seqs ({:d}) for {}".format(len(data["seqs"]), gs["TotEntries"], genome)
                )

            locTab = self.h5.create_table(
                "/Protein/Locus",
                genome,
                tablefmt.LocusTable,
                createparents=True,
                expectedrows=gs["TotEntries"] * 4,
            )

            cnt_missmatch_locus = 0
            cnt_genes = 0
            for nr in range(gs["TotEntries"]):
                eNr = data["off"] + nr + 1
                protTab.row["EntryNr"] = eNr
                protTab.row["OmaGroup"] = data["ogs"][nr]

                seqOff += self._add_sequence(data["seqs"][nr], protTab.row, seqArr, seqOff)
                cdnaOff += self._add_sequence(data["cdna"][nr], protTab.row, cdnaArr, cdnaOff, "CDNA")

                protTab.row["Chromosome"] = data["chrs"][nr]
                protTab.row["AltSpliceVariant"] = data["alts"][nr]
                protTab.row["OmaHOG"] = b" "  # will be assigned later
                protTab.row["CanonicalId"] = b" "  # will be assigned later
                if protTab.row["AltSpliceVariant"] == 0 or protTab.row["AltSpliceVariant"] == protTab.row["EntryNr"]:
                    cnt_genes += 1  # main isoforms of gene

                locus_str = data["locs"][nr]
                try:
                    locus_tab = loc_parser.parse(locus_str, eNr)
                    locTab.append(locus_tab)
                    len_cds = sum(z["End"] - z["Start"] + 1 for z in locus_tab)
                    if len_cds != protTab.row["CDNABufferLength"] - 1:
                        self.logger.debug(
                            "sum of exon lengths differ with cdna sequence for {}: {} vs {}".format(
                                eNr, len_cds, protTab.row["CDNABufferLength"] - 1
                            )
                        )
                        cnt_missmatch_locus += 1

                    protTab.row["LocusStart"] = locus_tab["Start"].min()
                    protTab.row["LocusEnd"] = locus_tab["End"].max()
                    protTab.row["LocusStrand"] = locus_tab[0]["Strand"]
                except ValueError as e:
                    self.logger.warning(e)
                protTab.row["SubGenome"] = data["subgenome"][nr].encode("ascii")
                protTab.row.append()
            protTab.flush()
            seqArr.flush()
            gs["TotGenes"] = cnt_genes
            gs.update()
            if cnt_missmatch_locus > 0:
                self.logger.warning("{} missmatches in exon-lengths compared to locus info".format(cnt_missmatch_locus))
            for n in (protTab, seqArr, locTab):
                if n.size_in_memory != 0:
                    self.logger.info(
                        "worte %s: compression ratio %3f%%" % (n._v_pathname, 100 * n.size_on_disk / n.size_in_memory)
                    )
        create_index_for_columns(protTab, "EntryNr", "MD5ProteinHash")
