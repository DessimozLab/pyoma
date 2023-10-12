from __future__ import division, print_function, unicode_literals

import collections
import functools
import io
import itertools
import json
import logging
import os
import re
import threading
import time
import warnings
from bisect import bisect_left
from builtins import chr, range, object, zip, bytes, str
from xml.etree import ElementTree as et
from typing import Union, Tuple, Optional, AnyStr, Set, Generator, Mapping

import dateutil
import fuzzyset
import networkx as nx
import numpy
import numpy.lib.recfunctions
import pandas
import pandas as pd
import pyopa
import tables
import tables.file as _tables_file
from Bio.UniProt import GOA
from datasketch import MinHash
from tqdm import tqdm

from .KmerEncoder import KmerEncoder
from .decorators import timethis, outdated_database_warning
from .exceptions import (
    InvalidTaxonId,
    DBVersionError,
    DBConsistencyError,
    DBOutdatedError,
    InvalidId,
    InvalidOmaId,
    UnknownIdType,
    UnknownSpecies,
    OutdatedHogId,
    Singleton,
    NoReprEntry,
    AmbiguousID,
    TooUnspecificQuery,
)
from .geneontology import GeneOntology, OntologyParser, GOAspect, FreqAwareGeneOntology
from .hoghelper import compare_levels, are_orthologous
from .hogprofile import Profiler
from .models import LazyProperty, KeyWrapper, ProteinEntry, Genome, HOG
from .suffixsearch import SuffixSearcher, SuffixIndexError
from .idmapper import (
    XRefNoApproximateIdMapper,
    XrefIdMapper,
    LinkoutIdMapper,
    GeneNameOrSymbolIdMapper,
    GeneNameAndMainDbIdMapper,
    UniProtIdMapper,
    DomainNameIdMapper,
)
from .. import version

logger = logging.getLogger(__name__)
warnings.filterwarnings(
    "ignore",
    category=tables.PerformanceWarning,
    message=r".*maximum recommended rowsize.*",
)

# Raise stack limit for PyOPA ~400MB
threading.stack_size(4096 * 100000)

# Global initialisations
GAF_VERSION = "2.1"


def count_elements(iterable):
    """return the number of elements in an iterator in the most efficient way.

    Be aware that for unbound iterators, this method won't terminate!
    :param iterable: an iterable object.
    """
    counter = itertools.count()
    collections.deque(zip(iterable, counter), maxlen=0)  # (consume at C speed)
    return next(counter)


_first_cap_re = re.compile("(.)([A-Z][a-z]+)")
_all_cap_re = re.compile("([a-z0-9])([A-Z])")


def to_snail_case(name):
    """function to convert from CamelCase to snail_case"""
    s1 = _first_cap_re.sub(r"\1_\2", name)
    return _all_cap_re.sub(r"\1_\2", s1).lower()


###
# here we monkeypatch pytables to be more thread friendly
# see http://www.pytables.org/cookbook/threading.html
class ThreadsafeFileRegistry(_tables_file._FileRegistry):
    lock = threading.RLock()

    @property
    def handlers(self):
        return self._handlers.copy()

    def add(self, handler):
        with self.lock:
            return super(ThreadsafeFileRegistry, self).add(handler)

    def remove(self, handler):
        with self.lock:
            return super(ThreadsafeFileRegistry, self).remove(handler)

    def close_all(self):
        with self.lock:
            return super(ThreadsafeFileRegistry, self).close_all()


class ThreadsafeFile(_tables_file.File):
    def __init__(self, *args, **kargs):
        with ThreadsafeFileRegistry.lock:
            super(ThreadsafeFile, self).__init__(*args, **kargs)

    def close(self):
        with ThreadsafeFileRegistry.lock:
            super(ThreadsafeFile, self).close()


@functools.wraps(tables.open_file)
def synchronized_open_file(*args, **kwargs):
    with ThreadsafeFileRegistry.lock:
        return _tables_file._original_open_file(*args, **kwargs)


# monkey patch the tables package
_tables_file._original_open_file = _tables_file.open_file
_tables_file.open_file = synchronized_open_file
tables.open_file = synchronized_open_file
_tables_file._original_File = _tables_file.File
_tables_file.File = ThreadsafeFile
tables.File = ThreadsafeFile
_tables_file._open_files = ThreadsafeFileRegistry()


# end monkey patch pytables
####


def count_rows_of_index_column_with_value(
    tab: Union[str, tables.Table], col: str, start, end=None
):
    idx = tab.colindexes[col]
    stop = end if end is not None else start
    return idx.search((start, stop))


def _get_limited_synteny_graph(
    edges_data: Generator[Tuple[int, int, dict], None, None],
    hog_tab: Union[tables.Table, numpy.ndarray, dict],
    ref_hog_idx: Optional[int],
    steps: Optional[int],
):
    G = nx.Graph()
    G.add_edges_from(edges_data)
    if ref_hog_idx is not None:
        try:
            neighbors = [ref_hog_idx] + [
                v for u, v in nx.bfs_edges(G, source=ref_hog_idx, depth_limit=steps)
            ]
        except nx.exception.NetworkXError as e:
            logger.error(
                "Trying to access a non-stored HOG in the ancestral synteny graph: "
                "hog_row: {}, Size of Graph: N={}, E={}".format(
                    ref_hog_idx, len(G), len(G.edges)
                )
            )
            raise DBConsistencyError(
                "Cannot find HogIdx {}  in Ancestral synteny graph".format(ref_hog_idx)
            )
        S = G.subgraph(neighbors)
    else:
        if isinstance(hog_tab, dict):
            G.add_nodes_from(hog_tab.keys())
        else:
            G.add_nodes_from(range(len(hog_tab)))
        S = G
    for n in S.nodes:
        hog = hog_tab[n]
        S.nodes[n].update(
            {
                "nr_members": int(hog["NrMemberGenes"]),
                "completeness_score": float(hog["CompletenessScore"]),
            }
        )
    return nx.relabel_nodes(S, lambda x: hog_tab[x]["ID"].decode())


def read_table_where(
    tab: tables.Table, query: str, condvars: Optional[Mapping] = None, field: str = None
) -> numpy.ndarray:
    """Reads data of a pytables table fulfilling the query as a numpy array.

    This is a more efficient implementation of the :py:method:`tables.Table.read_where`
    method.

    :param tab: The table to be read from
    :param query: the query
    :param condvars: mapping of condition variables of the query
    :param field: which column to be read.
                  If set to None, all columns are read, if set to 'nrow',
                  the row number is returned
    :return the data as a numpy array
    """
    it = tab.where(query, condvars=condvars)
    try:
        first = next(it)
    except StopIteration:
        return numpy.array([], dtype=tab.dtype)
    if field is None:
        fun = lambda row: row.fetch_all_fields()
        dtyp = tab.dtype
    elif field.lower == "nrow":
        fun = lambda row: row.nrow
        dtyp = numpy.int64
    else:
        fun = lambda row: row[field]
        dtyp = tab.dtype[field]

    data = numpy.fromiter(
        map(fun, itertools.chain([first], it)),
        dtype=dtyp,
    )
    return data


class Database(object):
    """This is the main interface to the oma database. Queries
    will typically be issued by methods of this object. Typically,
    the result of queries will be :py:class:`numpy.recarray` objects."""

    EXPECTED_DB_SCHEMA = "3.6"

    def __init__(self, db):
        if isinstance(db, str):
            logger.info("opening {} for read-only".format(db))
            self.db = tables.open_file(db, "r")
            self._close_fh = True
        elif isinstance(db, (tables.File, tables.file.File)):
            self.db = db
            self._close_fh = False
        else:
            raise ValueError(str(db) + " is not a valid database type")

        try:
            lib_version = self.db.get_node_attr("/", "pyoma_version")
        except AttributeError:
            lib_version = "<=0.7.0"
        logger.info(
            "pyoma library version: {}; db written with {}".format(
                version(), lib_version
            )
        )

        try:
            db_version = self.db.get_node_attr("/", "db_schema_version")
        except AttributeError:
            db_version = "1.0"

        logger.info("database version: {}".format(db_version))
        logger.info(
            "release: {}; release_char in HOG-IDs: '{}'".format(
                self.get_release_name(), self.release_char
            )
        )
        if db_version != self.EXPECTED_DB_SCHEMA:
            exp_tup = self.EXPECTED_DB_SCHEMA.split(".")
            db_tup = db_version.split(".")
            if db_tup[0] != exp_tup[0]:
                raise DBVersionError(
                    "Unsupported database version: {} != {} ({})".format(
                        db_version, self.EXPECTED_DB_SCHEMA, self.db.filename
                    )
                )
            else:
                logger.warning(
                    "outdated database version, but only minor version change: "
                    "{} != {}. Some functions might fail".format(
                        db_version, self.EXPECTED_DB_SCHEMA
                    )
                )
        self.db_schema_version = tuple(int(z) for z in db_version.split("."))
        self._on_close_notify = []

        try:
            self.seq_search = SequenceSearch(self)
        except DBConsistencyError as e:
            logger.exception(
                "Cannot load SequenceSearch. Any future call to seq_search will fail!"
            )
            self.seq_search = None
        self.id_resolver = IDResolver(self)
        self.id_mapper = IdMapperFactory(self)
        genomes = [Genome(self, g) for g in self.db.root.Genome.read()]
        self.tax = Taxonomy(
            self.db.root.Taxonomy.read(), genomes={g.ncbi_taxon_id: g for g in genomes}
        )
        try:
            self.desc_searcher = DescriptionSearcher(self)
        except SuffixIndexError:
            self.desc_searcher = None
        self.hog_profiler = None
        self._re_fam = None
        self.format_hogid = None
        self.mds = None
        self._set_hogid_schema()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def register_on_close(self, callback):
        self._on_close_notify.append(callback)
        logger.debug("registered call to {} on close".format(callback))

    def unregister_on_close(self, callback):
        try:
            self._on_close_notify.remove(callback)
            logger.debug("un-register {} on close".format(callback))
        except ValueError:
            logger.debug("could not un-register {}".format(callback))

    def close(self):
        for listener in self._on_close_notify:
            listener()
        if self._close_fh:
            self.get_hdf5_handle().close()

    @LazyProperty
    def gene_ontology(self):
        """returns GeneOntology object containing hierarchy
        of terms using the is_a and part_of relations. See
        :meth:`load_gene_ontology` to parametrize the
        creation of GeneOntology object."""
        return self.load_gene_ontology(GeneOntology)

    @LazyProperty
    def freq_aware_gene_ontology(self):
        gof = self.load_gene_ontology(factory=FreqAwareGeneOntology)
        try:
            go_count_tab = self.db.get_node("/Ontologies/GeneOntologyTermCounts")
            gof.cnts = {int(row["TermNr"]): int(row["Counts"]) for row in go_count_tab}
            gof.tot_cnts = [
                go_count_tab.attrs[z]
                for z in (
                    "total_molecular_function",
                    "total_biological_process",
                    "total_cellular_component",
                )
            ]
        except tables.NoSuchNodeError:
            pass
        return gof

    def load_gene_ontology(self, factory=None, rels=None):
        """Instantiate GeneOntology object

        By default, a GeneOntology object is returned based on
        the default relations (which are defined in :mod:`.gene_ontology`)

        The factory parameter allows to specify an subtype of
        GeneOntology, e.g. :class:`.gene_ontology.FreqAwareGeneOntology`,

        The rels parameter should be a list of relation strings that
        should be used as parents relations.

        :param factory: GeneOntology factory
        :param rels: list of rels for parent relations

        :returns: GeneOntology object"""
        try:
            fp = io.StringIO(
                self.db.root.Ontologies.GO.read().tobytes().decode("utf-8")
            )
        except tables.NoSuchNodeError:
            p = os.path.join(os.path.dirname(self.db.filename), "go-basic.obo")
            fp = open(p, "rt")
        if factory is None:
            factory = GeneOntology
        go = factory(OntologyParser(fp), rels=rels)
        go.parse()
        fp.close()
        return go

    def get_hdf5_handle(self):
        """return the handle to the database hdf5 file"""
        return self.db

    def get_conversion_date(self):
        """return the conversion end date from the DB attributes"""
        return dateutil.parser.parse(self.db.root._v_attrs["conversion_end"])

    def ensure_entry(self, entry):
        """This method allows to use an entry or an entry_nr.

        If necessary it will load the entry from the entry_nr,
        otherwise returning the same object again.

        :param entry: the entry_nr of a protein to be loaded or a
            protein entry."""
        try:
            t = entry["AltSpliceVariant"]
            return entry
        except (TypeError, AttributeError, IndexError):
            if isinstance(entry, (int, numpy.number)):
                return self.entry_by_entry_nr(entry)
            raise TypeError("Invalid type to retrieve an Entry")
        except Exception:
            raise TypeError("Invalid type to retrieve an Entry")

    def entry_by_entry_nr(self, entry_nr):
        """Returns the entry from the /Protein/Entries table
        corresponding to entry_nr.

        :param int entry_nr: a numeric identifier for the protein
            entry"""
        entry = self.db.root.Protein.Entries[entry_nr - 1]
        if entry["EntryNr"] != entry_nr:
            logger.warning(
                "EntryNr {} not at position {}. Using index instead".format(
                    entry_nr, entry_nr - 1
                )
            )
            entry = self.db.root.Protein.Entries.read_where(
                "EntryNr == {:d}".format(entry_nr)
            )
            if len(entry) != 1:
                raise ValueError(
                    "there are {} entries with entry_nr {}".format(len(entry), entry_nr)
                )
            entry = entry[0]
        return entry

    def _set_hogid_schema(self):
        """Determines the used HOG ID schema

        Some versions of the database have HOG IDs of the form
        "HOG:0000001" and others without the prefix (e.g. standalone)
        or with the prefix, but without padding. This method checks
        which schema is used and sets the appropriate member vars
        """
        re_id = re.compile(rb"(?P<prefix>HOG:)(?P<rel>[A-Z]+)?(?P<nr>\d+)")
        for entry in self.db.root.Protein.Entries:
            m = re_id.match(entry["OmaHOG"])
            if m is None:
                continue
            nr = m.group("nr")
            if len(nr) >= 7 and not nr.startswith(b"0"):
                continue  # a case where we cannot determine if padded nr
            is_padded = nr.startswith(b"0")
            prefix = m.group("prefix").decode()
            if prefix is None:
                prefix = ""
            rel = m.group("rel").decode() if m.group("rel") else ""
            if rel != self.release_char:
                raise DBConsistencyError(
                    "release_char does not match HOG-ids: {} vs {}".format(
                        self.release_char, rel
                    )
                )
            fmt = "{}{}{{:{}d}}".format(prefix, rel, "07" if is_padded else "")
            self._re_fam = re.compile(
                r"{}(?P<rel>[A-Z]+)?(?P<fam>\d+)".format(prefix).encode("ascii")
            )
            self.format_hogid = lambda fam: fmt.format(fam)
            logger.info(
                "setting HOG ID schema: re_fam: {}, hog_fmt: {}".format(
                    self._re_fam, fmt
                )
            )
            return
        raise DBConsistencyError("no protein in a hog")

    def all_proteins_of_genome(self, genome):
        """return all protein entries of a genome"""
        rng = self.id_mapper["OMA"].genome_range(genome)
        prot_tab = self.get_hdf5_handle().get_node("/Protein/Entries")
        return read_table_where(
            prot_tab, "(EntryNr >= {}) & (EntryNr <= {})".format(rng[0], rng[1])
        )

    def main_isoforms(self, genome):
        """returns the proteins that are the main isoforms of a genome.

        The main isoform is in this context the isoform that we used in OMA to
        infer the orthologs. It is the one variant that has the most alignment
        matches to all other gnomes.

        The `genome` parameter should be the UniProtSpeciesCode of the species of
        interest. If it is a numeric value, the genome parameter is interpreted
        as the protein entrynr. The method returns then the main isoforms for
        the species to which this protein belongs.

        :Note: OMA only predicts orthologs for the main isoform, so there is no
            difference if you work with only the main isoforms or all proteins of
            a genome in terms of orthologs.

        :param genome: UniProtSpeciesCode of the genome of interest, or a gene
                       number (EntryNr) from the genome of interest.
        """
        rng = self.id_mapper["OMA"].genome_range(genome)
        it = self._main_isoform_iter(rng[0], rng[1])
        try:
            first = next(it)
        except StopIteration:
            return numpy.array(
                [], dtype=self.get_hdf5_handle().get_node("/Protein/Entries").dtype
            )
        return numpy.fromiter(
            map(lambda x: x.fetch_all_fields(), itertools.chain([first], it)),
            dtype=first.table.dtype,
        )

    def _main_isoform_iter(self, entry_nr_low=None, entry_nr_high=None):
        query = []
        if entry_nr_low is not None:
            query.append("(EntryNr >= {})".format(entry_nr_low))
        if entry_nr_high is not None:
            query.append("(EntryNr <= {})".format(entry_nr_high))
        query.append("((AltSpliceVariant == EntryNr) | (AltSpliceVariant == 0))")
        query = " & ".join(query)
        prot_tab = self.get_hdf5_handle().get_node("/Protein/Entries")
        res_iter = prot_tab.where(query)
        return res_iter

    def count_main_isoforms(self, genome=None):
        """returns the number of main isoforms (i.e. the number
        of genes).

        If `genome` parameter is not `None`, the count is limited to the specified
        genome, otherwise all genomes in the database are used.

        The `genome` parameter should be the UniProtSpeciesCode of the species of
        interest. If it is a numeric value, the genome parameter is interpreted
        as the protein entrynr. The method returns then the main isoforms for
        the species to which this protein belongs.

        :param genome: UniProtSpeciesCode of the genome of interest, or a gene
                       number (EntryNr) from the genome of interest.
        """
        low, high = None, None
        if genome is not None:
            low, high = self.id_mapper["OMA"].genome_range(genome)
        return count_elements(self._main_isoform_iter(low, high))

    def get_splicing_variants(self, entry):
        e = self.ensure_entry(entry)
        if e["AltSpliceVariant"] == 0:
            return numpy.array([e], dtype=e.dtype)
        # TODO: create index on AltSpliceVariant column?!
        return read_table_where(
            self.get_hdf5_handle().get_node("/Protein/Entries"),
            "(EntryNr >= {:d}) & (EntryNr < {:d}) & (AltSpliceVariant == {:d})".format(
                e["EntryNr"] - 100, e["EntryNr"] + 100, e["AltSpliceVariant"]
            ),
        )

    def get_extant_synteny_graph(
        self,
        genome,
        center_entry: Optional[Union[int, str]] = None,
        window: Optional[int] = None,
    ):
        all_genes = self.main_isoforms(genome)
        all_genes.sort(order=["Chromosome", "LocusStart"])
        gs = self.id_mapper["OMA"].genome_of_entry_nr(all_genes[0]["EntryNr"])
        if center_entry is not None:
            if window is None:
                window = 7
            if isinstance(center_entry, str) and center_entry.startswith("HOG:"):
                k = numpy.where(all_genes["OmaHOG"] == center_entry.encode("utf-8"))[0]
                if len(k) == 1:
                    en = all_genes[k]
                else:
                    raise InvalidId(
                        f"{center_entry} unknown for {gs['UniProtSpeciesCode']}"
                    )
            else:
                en = self.ensure_entry(self.id_resolver.resolve(center_entry))
                if (
                    en["AltSpliceVariant"] > 0
                    and en["AltSpliceVariant"] != en["EntryNr"]
                ):
                    en = self.ensure_entry(en["AltSpliceVariant"])
            try:
                idx = numpy.where(all_genes["EntryNr"] == en["EntryNr"])[0][0]
                ref = all_genes[idx]
            except IndexError:
                identified_genome = gs["UniProtSpeciesCode"].decode()
                raise InvalidId(
                    f"Center gene {center_entry} is not a valid for genome {identified_genome}"
                )
            all_genes = all_genes[idx - window : idx + window + 1]
            all_genes = all_genes[
                numpy.where(all_genes["Chromosome"] == ref["Chromosome"])
            ]
        oma_ids = list(
            map(
                lambda x: f"{gs['UniProtSpeciesCode'].decode()}{x - gs['EntryOff']:05d}",
                all_genes["EntryNr"],
            )
        )

        def node_data_generator(genes):
            for id_, g in zip(oma_ids, genes):
                yield (
                    id_,
                    {
                        "chromosome": g["Chromosome"].decode(),
                        "start": int(g["LocusStart"]),
                        "strand": "+" if g["LocusStrand"] > 0 else "-",
                        "hog_id": g["OmaHOG"].decode(),
                    },
                )

        G = nx.Graph()
        G.add_nodes_from(node_data_generator(all_genes))
        for i in range(len(all_genes) - 1):
            if all_genes[i]["Chromosome"] == all_genes[i + 1]["Chromosome"]:
                G.add_edge(oma_ids[i], oma_ids[i + 1], weight=1)
        return G

    def _get_vptab(self, entry_nr):
        return self._get_pw_tab(entry_nr, "VPairs")

    def _get_pw_tab(self, entry_nr, subtab):
        genome = (
            self.id_mapper["OMA"]
            .genome_of_entry_nr(entry_nr)["UniProtSpeciesCode"]
            .decode()
        )
        return self.db.get_node("/PairwiseRelation/{}/{}".format(genome, subtab))

    def nr_ortholog_relations(self, entry_nr):
        """returns the number of orthologous relations for a given entry.

        The return value is a numpy.void tuple with all the gene centric orthology/paralogy
        counts. The named fields are:

          - *NrPairwiseOrthologs*: nr of pairwise orthologs as computed by OMA
          - *NrHogInducedPWOrthologs*: nr of induced pairwise orthologs according
            HOG inference process
          - *NrHogInducedPWParalogs*: nr of induced pairwise paralogs
          - *NrOMAGroupOrthologs*: nr of orthologs according to the OMA Groups
          - *NrAnyOrthologs*: sum of all orthologs according to pairwise orthologs,
            induced pairwise orthologs and OMA Group based orthologs.

        :param int entry_nr: entry number to check
        """
        tab = self.db.get_node("/Protein/OrthologsCountCache")
        return tab[entry_nr - 1]

    def count_vpairs(self, entry_nr):
        """number of pairwise orthologs

        :see_also: :meth:`nr_ortholog_relations`"""
        try:
            return int(self.nr_ortholog_relations(entry_nr)["NrPairwiseOrthologs"])
        except tables.NoSuchNodeError:
            # fallback in case Cache does not exist yet
            vptab = self._get_vptab(entry_nr)
            try:
                cnt = count_elements(vptab.where("(EntryNr1=={:d})".format(entry_nr)))
            except (TypeError, ValueError):
                cnt = 0
            return cnt

    def count_homoeologs(self, entry_nr):
        """number of homoeologs of a given protein entry"""
        pwtab = self._get_pw_tab(entry_nr, "within")
        homolog_typ_nr = pwtab.get_enum("RelType")["homeolog"]
        try:
            cnt = count_elements(
                pwtab.where(
                    "(EntryNr1=={:d}) & (RelType == {:d})".format(
                        entry_nr, homolog_typ_nr
                    )
                )
            )
        except (TypeError, ValueError):
            cnt = 0
        return cnt

    def _get_pw_data(
        self, entry_nr, tab, target_range=None, typ_filter=None, extra_cols=None
    ):
        if isinstance(entry_nr, tuple):
            query = "(EntryNr1 >= {:d}) & (EntryNr1 <= {:d})".format(*entry_nr)
        else:
            query = "(EntryNr1 == {:d})".format(entry_nr)
        if target_range is not None:
            query += " & (EntryNr2 >= {:d}) & (EntryNr2 <= {:d})".format(*target_range)
        if typ_filter is not None:
            query += " & (RelType == {:d})".format(typ_filter)
        dat = read_table_where(tab, query)
        typ = tab.get_enum("RelType")
        cols = ["EntryNr1", "EntryNr2", "Score", "Distance"]
        if extra_cols is not None:
            cols.extend(extra_cols)
        res = numpy.lib.recfunctions.append_fields(
            dat[cols],
            names="RelType",
            data=[typ(x) for x in dat["RelType"]],
            usemask=False,
        )
        return res

    def get_vpairs(self, entry_nr):
        """returns the verified pairs of a query protein.

        This method returns an instance of a :class:`numpy.recarray` class
        containing the verified pairs of a query protein entry.
        The returned array contains columns with EntryNr1 and EntryNr2 to
        identify the pair together with RelType (indicating the subtype of
        orthology), the alignment score and the distance. The score and
        distance will be set to -1 if unknown.

        :param int entry_nr: the numeric entry_nr of the query protein."""
        vp_tab = self._get_vptab(entry_nr)
        return self._get_pw_data(entry_nr, vp_tab)

    def get_vpairs_between_species_pair(self, genome1, genome2):
        """Returns all the pairwise orthologs between two species.

        This method returns a numpy array of the same dtype as
        :meth:`get_vpairs` does. As input, two genomes present in the
        OMA database need to be passed, either as rows from from the
        GenomeTable or with the UniProt mnemonic species code.
        The pairwise orthologs are symetric, so the order of genome1
        and genome2 does not have any impact (except of the order of
        EntryNr1 and EntryNr2)

        :param genome1: first genome
        :type genome1: :class:`numpy.void`, :class:`pyoma.browser.models.Genome`, str
        :param genome2: second genome
        :type genome2: :class:`numpy.void`, :class:`pyoma.browser.models.Genome`, str
        :returns array with pairwise orthologs between genome1 and genome2
        :rtype: `numpy.ndarray` of dtype :class:`.tablefmt.PairwiseRelationTable`

        """

        def to_genome(g) -> Genome:
            if isinstance(g, Genome):
                return g
            else:
                return Genome(self, g)

        g1, g2 = map(to_genome, (genome1, genome2))
        vptab = self._get_vptab(g1.entry_nr_offset + 1)
        return self._get_pw_data(
            entry_nr=(g1.entry_nr_offset + 1, g1.entry_nr_offset + g1.nr_entries),
            target_range=(g2.entry_nr_offset + 1, g2.entry_nr_offset + g2.nr_entries),
            tab=vptab,
        )

    def get_within_species_paralogs(self, entry_nr):
        """returns the within species paralogs of a given entry

        This method returns a :class:`numpy.recarray` instance
        containing the close paralogs. Close paralogs are within
        species paralogs that are inparalogs to at least one
        ortholog of the query gene in OMA.

        The returned array contains columns with EntryNr1 and EntryNr2 to
        identify the pair together with RelType (indicating the subtype of
        paralogy), the alignment score and the distance. The score and
        distance will be set to -1 if unknown.

        :param int entry_nr: the numeric entry_id of the query protein"""
        within_species_paralogs = self._get_pw_tab(entry_nr, "within")
        return self._get_pw_data(entry_nr, within_species_paralogs)

    def get_homoeologs(self, entry_nr):
        within_species = self._get_pw_tab(entry_nr, "within")
        homolog_typ_nr = within_species.get_enum("RelType")["homeolog"]
        return self._get_pw_data(
            entry_nr,
            within_species,
            typ_filter=homolog_typ_nr,
            extra_cols=["SyntenyConservationLocal", "Confidence"],
        )

    def get_hog_induced_pairwise_orthologs(self, entry):
        """This method retrieves the hog induced pairwise orthologs

        The induced pairwise orthologs are in general not equal to
        the vpairs. The method returns a numpy array with the entries
        that are orthologous to the query entry

        :param entry: entry or entry_nr of the query protein"""
        entry = self.ensure_entry(entry)
        genome_entry_range = self.id_mapper["OMA"].genome_range(entry["EntryNr"])

        def is_orthologous(a, b):
            """genes are orthologs if their HOG id have a common prefix that is
            either the base id of the family or the prefix does not end with
            a subfamily number, ie. not a digit as common prefix. See LOFT paper
            for details on encoding."""
            if a["EntryNr"] == b["EntryNr"]:
                return False
            prefix = os.path.commonprefix((a["OmaHOG"], b["OmaHOG"])).decode()
            if "." in prefix and prefix[-1].isdigit():
                return False
            # count number of genes in query genome that are co-orthologs (== having the prefix)
            cnts = numpy.char.startswith(
                hogids_of_genes_in_query_genome, prefix.encode("utf-8")
            ).sum()
            return cnts

        try:
            fam = self.hog_family(entry)
            hog_member = self.member_of_fam(fam)
        except Singleton:
            # an empty fetch
            hog_member = self.db.root.Protein.Entries[0:0]
        hogids_of_genes_in_query_genome = hog_member[
            numpy.where(
                (hog_member["EntryNr"] >= genome_entry_range[0])
                & (hog_member["EntryNr"] < genome_entry_range[1])
            )
        ]["OmaHOG"]
        query_genome_genes_cnt = numpy.array(
            [is_orthologous(entry, hog_member[i]) for i in range(len(hog_member))],
            dtype="i4",
        )
        mask = numpy.asarray(query_genome_genes_cnt, bool)
        target_genomes = [
            self.id_mapper["OMA"].genome_of_entry_nr(o["EntryNr"])["NCBITaxonId"]
            for o in hog_member[mask]
        ]
        targ_genome_genes_cnt = collections.Counter(target_genomes)
        induced_orthologs = numpy.lib.recfunctions.append_fields(
            hog_member[mask],
            names="RelType",
            data=[
                "{}:{}".format(
                    "m" if cnt_a > 1 else 1, "n" if targ_genome_genes_cnt[b] > 1 else 1
                ).encode("utf-8")
                for cnt_a, b in zip(query_genome_genes_cnt[mask], target_genomes)
            ],
        )
        return induced_orthologs

    def get_hog_induced_pairwise_paralogs(self, entry, start=0, stop=None):
        """This method retrieves the hog induced pairwise paralogs

        The method returns a numpy array with the entries
        that are paralogous to the query entry. In addition to the normal
        ProteinEntry, the array contains an additional column with the
        taxonomic level when the duplication occurred.

        :param entry: entry or entry_nr of the query protein"""
        entry = self.ensure_entry(entry)
        genome = self.id_mapper["OMA"].genome_of_entry_nr(entry["EntryNr"])
        lineage = self.tax.get_parent_taxa(genome["NCBITaxonId"])["Name"]
        lineage = numpy.fromiter(
            filter(lambda x: x in self.tax.all_hog_levels, lineage), dtype=lineage.dtype
        )
        lineage_sorter = numpy.argsort(lineage)
        try:
            fam = self.hog_family(entry)
        except Singleton:
            # an empty fetch
            hog_member = self.db.root.Protein.Entries[0:0]
        else:
            levels = self.filter_fam_from_hoglevel(fam)
            # only keep the levels that are on the lineage to the query genome
            levels = levels[numpy.isin(levels["Level"], lineage)]
            hog_member = self._members_of_hog_id(self.format_hogid(fam))

        def is_paralogous(a, b):
            """genes are orthologs if their HOG id have a common prefix that is
            either the base id of the family or the prefix does not end with
            a subfamily number, ie. not a digit as common prefix. See LOFT paper
            for details on encoding."""
            if a["EntryNr"] == b["EntryNr"]:
                return False
            prefix = os.path.commonprefix((a["OmaHOG"], b["OmaHOG"])).decode()
            if "." in prefix and prefix[-1].isdigit():
                # gene is paralog. find MRCA in taxonomy of common HOGid prefix
                k = prefix.rfind(".")
                hog_id = prefix[:k].encode("ascii")
                cand_levels = levels[numpy.where(levels["ID"] == hog_id)]
                sortidx = lineage.searchsorted(
                    cand_levels["Level"], sorter=lineage_sorter
                )
                lin_idx = numpy.take(lineage_sorter, sortidx, mode="clip")
                mask = lineage[lin_idx] == cand_levels["Level"]
                # we take the first diverged lineage, meaning the duplication happened
                # on the branch to that level.
                return lineage[lin_idx[mask].min() - 1]
            return None

        def filter_candidates(entry, candidates):
            for cand in candidates:
                lev = is_paralogous(entry, cand)
                if lev:
                    yield tuple(cand) + (lev,)

        ext_prot_dtype = numpy.dtype(entry.dtype.descr + [("DivergenceLevel", "S255")])
        paralogs = numpy.fromiter(
            itertools.islice(filter_candidates(entry, hog_member), start, stop),
            dtype=ext_prot_dtype,
        )
        return paralogs

    def get_hog_induced_orthologs_between_genome_pair(
        self, g1, g2, return_cardinality=False, return_omagroup=False
    ):
        """Returns the HOG-induced pairwise orthologs between two genomes.

        The method returns a numpy structured array with EntryNr1, EntryNr2 and HogID columns.
        If `return_cardinaltiy` is set to True, a `RelType` column is added which returns the
        cardinality of the relation (1:1, 1:n, m:1, m:1).
        If `return_omagroup` is set to True, the returned array also has a column `OmaGroup` the
        contains the numeric oma group if the two proteins belong to the same oma group. Otherwise,
        the value is set to NaN.

        :param g1: first genome
        :type g1: :class:`numpy.void`, :class:`pyoma.browser.models.Genome`, str
        :param g2: second genome
        :type g2: :class:`numpy.void`, :class:`pyoma.browser.models.Genome`, str
        :param return_cardinality: whether or not to return the cardinality, i.e. 1:1, 1:n, m:1, m:n
        :type return_cardinality: bool
        :param return_omagroup: whether or not to return the oma group nr in case both proteins
                                belong to the same oma group. If not, NaNs will be returned.
        :type return_omagroup: bool
        :returns: array with hog induced pairwise orthologs between genome1 and genome2
        """

        def load_prot_dataframe(g):
            genome = g if isinstance(g, Genome) else Genome(self, g)
            prots = pandas.DataFrame(
                numpy.sort(
                    self.main_isoforms(genome.uniprot_species_code), order="OmaHOG"
                )
            )
            col_keep = {"EntryNr", "OmaHOG", "OmaGroup", "CanonicalId"}
            prots.drop(columns=set(prots.columns) - col_keep, inplace=True)
            prots["Fam"] = prots["OmaHOG"].apply(
                lambda x: 0 if x == b"" else self.parse_hog_id(x)
            )
            prots.drop(
                prots[prots.Fam == 0].index, inplace=True
            )  # remove prots that don't belong to a hog
            return prots

        prots1 = load_prot_dataframe(g1)
        prots2 = load_prot_dataframe(g2)
        joined_fam = prots1.join(
            prots2.set_index("Fam"), on="Fam", how="inner", lsuffix="1", rsuffix="2"
        )
        joined_fam["HogID"] = joined_fam.apply(
            lambda row: are_orthologous(
                row["OmaHOG1"], row["OmaHOG2"], return_group=True
            ),
            axis=1,
        )
        orthologs = joined_fam[joined_fam["HogID"] != False].copy()
        orthologs["HogID"] = orthologs["HogID"].str.encode("ascii")
        return_cols = ["EntryNr1", "EntryNr2", "HogID"]
        if return_cardinality:
            e1_cnts = orthologs["EntryNr1"].value_counts()
            e2_cnts = orthologs["EntryNr2"].value_counts()
            orthologs.loc[:, "RelType"] = orthologs.apply(
                lambda r: (
                    ("1" if e2_cnts[r.EntryNr2] == 1 else "m")
                    + ":"
                    + ("1" if e1_cnts[r.EntryNr1] == 1 else "n")
                ).encode("ascii"),
                axis=1,
            )
            return_cols.append("RelType")
        if return_omagroup:
            orthologs.loc[:, "OmaGroup"] = orthologs.apply(
                lambda r: r.OmaGroup1 if r.OmaGroup1 == r.OmaGroup2 else numpy.nan,
                axis=1,
            )
            return_cols.append("OmaGroup")
        orthologs = orthologs[return_cols]
        typ = orthologs.dtypes
        return numpy.array(
            [tuple(x) for x in orthologs.values], dtype=list(zip(typ.index, typ))
        )

    def neighbour_genes(self, entry_nr, window=1):
        """Returns neighbor genes around a query gene.

        This method returns a tuple containing a numpy recarray with
        gene entries located around the query gene, and an index
        pointing to the query gene. The genes are sorted according to
        their position on the chromosome.

        The *windows* parameter specifies the number of genes up- and
        downstream of the query gene that should be reported. Note
        that the actual number can be smaller if the query gene is close
        to a chromosome start or end.

        :param entry_nr: the entry number of the query gene
        :param window: the number of neighboring genes on each
                           side to return"""
        if window <= 0 or not isinstance(window, int):
            raise ValueError("windows parameters must be a positive integer value")

        dat = self.entry_by_entry_nr(entry_nr)
        target_chr = dat["Chromosome"]
        genome_range = self.id_mapper["OMA"].genome_range(entry_nr)
        f = 5
        condvars = {
            "enr_max": min(genome_range[1], entry_nr + f * window),
            "enr_min": max(genome_range[0], entry_nr - f * window),
            "chr": target_chr,
        }
        data = read_table_where(
            self.db.root.Protein.Entries,
            "(EntryNr >= enr_min) & (EntryNr <= enr_max) & "
            "(Chromosome == chr) & "
            "((AltSpliceVariant == 0) |"
            " (AltSpliceVariant == EntryNr))",
            condvars=condvars,
        )
        data.sort(order=["LocusStart"])
        idx = int((data["EntryNr"] == entry_nr).nonzero()[0])
        res = data[max(0, idx - window) : idx + window + 1]
        idx = int((res["EntryNr"] == entry_nr).nonzero()[0])
        return res, idx

    def parse_hog_id(self, hog_id):
        hog_id = hog_id if isinstance(hog_id, bytes) else hog_id.encode("ascii")
        m = self._re_fam.match(hog_id)
        if m is not None:
            rel = m.group("rel")
            if rel is None:
                rel = b""
            if rel.decode() != self.release_char:
                raise OutdatedHogId(hog_id)
            return int(m.group("fam"))
        else:
            raise ValueError("invalid hog id format")

    def hog_family(self, entry):
        entry = self.ensure_entry(entry)
        m = self._re_fam.match(entry["OmaHOG"])
        if m is None:
            raise Singleton(entry)
        return int(m.group("fam"))

    def hog_levels_of_fam(self, fam_nr, deduplicate_and_decode=False):
        """get all taxonomic levels covered by a family.

        The family corresponds to the toplevel numeric id of a HOG,
        i.e. for HOG:002421 the fam_nr should be 2421. If a HOG
        covers a certain level more than once, it will be returned
        several times, unless `deduplicate_and_decode` is set to True.

        :param int fam_nr: the numeric id of the family (== Toplevel HOG)

        :param bool deduplicate_and_decode: decode the encoded levels and
            return only the unique levels as a frozenset(string).
            Added in version 0.8.0
        """
        levels = self.filter_fam_from_hoglevel(fam_nr, field="Level")
        if deduplicate_and_decode:
            levels = frozenset(x.decode() for x in frozenset(levels))
        return levels

    @functools.lru_cache(maxsize=128)
    def filter_fam_from_hoglevel(self, fam_nr, field=None):
        t0 = time.time()
        hoglevel_tab = self.db.get_node("/HogLevel")
        try:
            fam_idx = self.db.get_node("/HogLevel_fam_lookup")
            levels = hoglevel_tab.read(*fam_idx[fam_nr], field=field)
        except IndexError:
            # dummy read that returns empty list of same dtype
            levels = hoglevel_tab.read(0, 0, field=field)
        except tables.NoSuchNodeError:
            # fall back to index based search
            levels = self.db.root.HogLevel.read_where(
                "(Fam=={})".format(fam_nr), field=field
            )
        logger.debug(
            "retrieving levels for family {:d} took {:.7f} sec".format(
                fam_nr, time.time() - t0
            )
        )
        return levels

    def get_subhogs(
        self, hog_id, level=None, include_subids=False, include_leaf_levels=True
    ):
        """Get all the (sub)hogs for a given hog_id

        This method returns an generator of :class:`models.HOG` instances,
        that have a certain (sub)hog id. If `include_subids` is set to
        False (default), only the (sub)HOGs with exactly the query
        hog_id will be returned, i.e. a set of taxonomic ranges for
        which no duplication occurred in between for this HOG.

        The method returns a generator of :class:`models.HOG` instances.

        :param (bytes, str) hog_id: the hog_id of interest

        :param str level: the root level, defaults to the root of the subhog

        :param bool include_subids: whether to include suhogs
                                    that originated after a duplication
                                    in the query (sub)HOG. defaults to False

        :param bool include_leaf_levels: whether to include the level
                                    of the extant species. defaults to True.
                                    Leaf levels have by definition only one
                                    member and are thus not of limited interest
                                    in most situations

        :returns: generator of HOG instances
        :rtype: :class:`models.HOG`

        :see_also: :meth:`get_hog` that returns a single HOG instance
            for a specific level or the root level one for a specific HOG id.
        """
        if level is not None and level not in ("LUCA", b"LUCA"):
            try:
                subtax = self.tax.get_subtaxonomy_rooted_at(level, collapse=False)
                if not include_leaf_levels:
                    idx = subtax.tax_table["NCBITaxonId"].searchsorted(
                        subtax.tax_table["ParentTaxonId"], sorter=subtax.taxid_key
                    )
                    internal_node_idx = set(
                        pos
                        for i, pos in enumerate(idx)
                        if 0 <= pos < len(subtax.tax_table)
                        and subtax.tax_table["NCBITaxonId"][pos]
                        == subtax.tax_table["ParentTaxonId"][i]
                    )
                    internal_node_idx = numpy.fromiter(
                        internal_node_idx, dtype="i4", count=len(internal_node_idx)
                    )
                    subtax.tax_table = subtax.tax_table[internal_node_idx]
                children = set(n["Name"] for n in subtax.tax_table)
                children.remove(
                    level.encode("utf-8") if isinstance(level, str) else level
                )
            except KeyError as e:
                raise ValueError(f"invalid level: {level}") from e
        elif not include_leaf_levels:
            children = self.tax.all_hog_levels

        if include_subids:
            condvars = {
                v: k
                for v, k in zip(
                    ("hogid_low", "hogid_high"), self._hog_lex_range(hog_id)
                )
            }
            query = "(ID >= hogid_low) & (ID < hogid_high)"
        else:
            query = "ID == hogid"
            condvars = {"hogid": hog_id}

        for row in self.db.root.HogLevel.where(query, condvars=condvars):
            hog = row.fetch_all_fields()
            if (
                (level is None or level in ("LUCA", b"LUCA")) and include_leaf_levels
            ) or hog["Level"] in children:
                yield HOG(self, hog)

    def get_subhogids_at_level(self, fam_nr, level):
        """get all the hog ids within a given family at a given taxonomic
        level of interest.

        After a duplication in an ancestor lineage, there exists multiple
        sub-hogs for any taxonomic level after the duplication. This method
        allows to get the list of hogids at the requested taxonomic level.

        E.g. assume in family 1 (HOG:0000001) there has been a duplication
        between Eukaryota and Metazoa. this method would return for
        get_subhogids_at_level(1, 'Eukaryota') --> ['HOG:0000001']
        and for
        get_subhogids_at_level(1, 'Metazoa') --> ['HOG:0000001.1a', 'HOG:0000001.1b']

        :note::
        There is also the method :meth:`get_subhogs_at_level` which returns all
        information stored in the HogLevel table, not only the HOG id.

        :see_also:: :meth:`get_subhogs_at_level`

        :param fam_nr: the numeric family id
        :param level: the taxonomic level of interest"""
        return self.get_subhogs_at_level(fam_nr, level, field="ID")

    def get_subhogs_at_level(self, fam_nr, level, field=None):
        """get all the hogs within a given family at a given taxonomic
        level of interest.

        After a duplication in an ancestor lineage, there exists multiple
        sub-hogs for any taxonomic level after the duplication. This method
        allows to get the list of hogs at the requested taxonomic level.

        :note:
        There is also the method get_subhogids_at_level which returns only
        the hogids for the chosen level, not all information on the HOGs.

        :param fam_nr: the numeric family id
        :param level: the taxonomic level of interest
        :param field: name of the column to be returned. If not set, all columns
            will be returned."""
        lev = level if isinstance(level, bytes) else level.encode("ascii")
        return read_table_where(
            self.db.root.HogLevel,
            query="(Fam==fam) & (Level == lev)",
            condvars={"fam": fam_nr, "lev": lev},
            field=field,
        )

    def get_parent_hogs(self, hog_id, level=None):
        """return an array of parent (deeper) HOGs.

        The returned array is sorted from deepest to shallowest parent
        HOG of the query HOG. The query HOG can be either defined just
        by its ID, in which case the deepest level of that (sub)hog is
        assumed, or a valid hog_id/level combination.

        At the last index, the query HOG will be present, so unless
        an invalid HOG is selected, at least one HOG will be returned.

        The method returns a list of :class:`models.HOG` instances

        :param (str,bytes) hog_id: a valid hog_id
        :param str level: optional level of reference HOG"""
        if isinstance(hog_id, str):
            hog_id = hog_id.encode("ascii")
        ref_hog = self.get_hog(hog_id, level=level)
        parent_pos = {b"LUCA": -1}
        if ref_hog["Level"] != b"LUCA":
            taxnode_of_level = self.tax.get_taxnode_from_name_or_taxid(ref_hog["Level"])
            parent_taxnodes = self.tax.get_parent_taxa(
                taxnode_of_level[0]["NCBITaxonId"]
            )
            parent_pos.update(
                {lev: pos for pos, lev in enumerate(parent_taxnodes["Name"][::-1])}
            )
        query = "(Fam == fam) & (ID <= hogid)"
        condvars = {"fam": ref_hog["Fam"], "hogid": hog_id}
        parent_hogs = [] * len(parent_pos)
        for row in self.db.root.HogLevel.where(query, condvars=condvars):
            hog = row.fetch_all_fields()
            if not hog_id.startswith(hog["ID"]):
                continue
            if len(hog_id) > len(hog["ID"]) and hog_id[len(hog["ID"])] != ord("."):
                continue
            if hog["Level"] in parent_pos:
                parent_hogs.append(HOG(self, hog))
        parent_hogs.sort(key=lambda x: parent_pos[x._hog["Level"]])
        return parent_hogs

    def _members_of_hog_id(self, hog_id):
        hog_range = self._hog_lex_range(hog_id)
        it = self.db.root.Protein.Entries.where(
            "(low <= OmaHOG) & (OmaHOG < high)",
            condvars={k: v for k, v in zip(("low", "high"), hog_range)},
        )
        # we need to filter them in case there are many (>26) paralog clusters,
        # in which case the hog_id need to be followed by a '.'
        hog_id_len = len(hog_range[0])

        def is_same_subhog(row):
            return (
                len(row["OmaHOG"]) == hog_id_len
                or chr(row["OmaHOG"][hog_id_len]) == "."
            )

        for row in it:
            if is_same_subhog(row):
                yield row.fetch_all_fields()

    def member_of_hog_id(self, hog_id, level=None):
        """return an array of protein entries which belong to a given hog_id.

        E.g. if hog_id = 'HOG122.1a', the method returns all the proteins that
        have either exactly this hog id or an inparalogous id such a HOG122.1a.4b.2a

        If you are only interested in the members of a specific lineage (identified
        through its taxonomic range), you can pass the taxonomic range as an
        additional argument. Only the proteins of genomes belonging to this clade
        will be returned. Otherwise, all proteins with having this specific hog_id
        will be returned.

        :param str hog_id: the requested hog_id.
        :param level: the taxonomic level of interest
        :type level: str or None

        :return: a numpy.array with the protein entries belonging to the requested hog.
        :rtype: :class:`numpy.ndarray`

        :Note: Even if you obtained a certain hog_id using
               :py:meth:`get_subhogids_at_level`
               using a certain level, if you do not specify the level in
               :meth:`member_of_hog_id` again, you will likely get proteins from other
               clades. Only if it happens that the deepest level of the hog_id
               coincides with the taxonomic range of interest, the two will be identical.
        """
        it = self._members_of_hog_id(hog_id)
        try:
            first = next(it)
        except StopIteration:
            return numpy.array([], dtype=self.db.root.Protein.Entries.dtype)

        members = numpy.fromiter(itertools.chain([first], it), dtype=first.dtype)
        if level is not None and level != "LUCA":
            keep = numpy.array(
                [
                    level.encode("ascii")
                    in self.tax.get_parent_taxa(
                        self.id_mapper["OMA"].genome_of_entry_nr(enr)["NCBITaxonId"]
                    )["Name"]
                    for enr in members["EntryNr"]
                ],
                dtype=bool,
            )
            members = members[keep]
        return members

    def iter_members_of_hog_id(self, hog_id, start=0, stop=None, step=1):
        """iterates over all proteins that belong to a specific hog_id.

        A hog_id might be an ID of the following form: HOG:0000212.1a
        This method will yield all proteins in the form of
        :class:`ProteinEntry` instances that are part of this hog_id.

        The paramters start, stop and step are passed to the
        :py:func:`itertools.islice` function and can be used to access
        only a subset of the members, e.g. for pageination.

        :param str hog_id: the requested HOG ID.
        :return: :py:class:`ProteinEntry` objects
        :rtype: iter(:class:`ProteinEntry`)"""
        members = self._members_of_hog_id(hog_id)
        for row in itertools.islice(members, start, stop, step):
            yield ProteinEntry(self, row)

    def member_of_fam(self, fam):
        """returns an array of protein entries which belong to a given fam"""
        if not isinstance(fam, (int, numpy.number)):
            raise ValueError("expect a numeric family id")
        return self.member_of_hog_id(self.format_hogid(fam))

    def hog_members(self, entry, level):
        """get hog members with respect to a given taxonomic level.

        The method will return a list of protein entries that are all
        member of the same hog with respect to the taxonomic range
        of interest.

        :param entry: an entry or entry_nr of a query protein
        :param level: the taxonomic level of interest"""
        query = self.ensure_entry(entry)
        members = self.hog_members_from_hog_id(query["OmaHOG"], level)
        if query not in members:
            raise ValueError("Level '{0:s}' undefined for query gene".format(level))
        return members

    def hog_members_from_hog_id(self, hog_id, level):
        """get hog members with respect to a given taxonomic level.

        The method will return a list of protein entries that are all
        member of the same hog with respect to the taxonomic range
        of interest.

        :param bytes hog_id: the query hog id
        :param str level: the taxonomic level of interest"""
        try:
            hog = next(self.iter_hogs_at_level(hog_id, level))
        except StopIteration:
            raise ValueError('Level "{0:s}" undefined for query gene'.format(level))
        # get the entries which have this hogid (or a sub-hog)
        members = self.member_of_hog_id(hog["ID"], level=level)
        return members

    def iter_hogs_at_level(self, hog_id, level):
        """yields the (sub- or parent-) hogs for a given level.

        This method yields the hogs as numpy.array instances that
        match the requested level exactly and match the hog_id as
        closely as possible:

        For a more general hog_id, the method returns all the subhogs
        that exist at the requested level (similar to
        :meth:`get_subhogs_at_level`, but here respecting a potential
        subhog part).

        For a hog_id which is more specific than the hog_id at the
        requested level, only the parent hog as defined at the level
        is returned.

        :Note:
        If the hog_id is known to match the level, the method
        :meth:`get_hog` can be used as it is a bit faster.

        :see_also: :meth:`get_subhogs_at_level`, :meth:`get_hog`
        """
        if isinstance(hog_id, str):
            hog_id = hog_id.encode("ascii")

        condvars = {"fam": self.parse_hog_id(hog_id)}
        query = "(Fam == fam)"
        try:
            anc_node = self._ancestral_node(level)
            hog_tab = anc_node.Hogs
        except DBConsistencyError:
            hog_tab = self.db.root.HogLevel
            query += " & (Level == lev)"
            condvars["lev"] = level.encode("utf-8")

        for hog in hog_tab.where(query, condvars=condvars):
            if hog_id.startswith(hog["ID"]):
                if hog_id == hog["ID"] or chr(hog_id[len(hog["ID"])]) == ".":
                    yield hog.fetch_all_fields()
            elif hog["ID"].startswith(hog_id):
                yield hog.fetch_all_fields()

    def iter_hog_at_level(self, hog_id, level):
        """deprecated method. use iter_hogs_at_level"""
        yield from self.iter_hogs_at_level(hog_id=hog_id, level=level)

    def get_hog(self, hog_id, level=None, field=None, tab=None):
        """Retrieve the one relevant HOG for a certain hog-id.

        If a level is provided, returns the (sub)hog at this level, otherwise
        it will return the deepest (sub)hog for that ID. Note that for this
        method the hog_id must match the expected subhog exactly. That means
        the method does not try to identify the correct hog_id for the requested
        level. In case the hog_id needs to be potentially adjusted given a level
        of interest, you should use :meth:`iter_hogs_at_level`.

        :see_also: :meth:`iter_hogs_at_level`, :meth:`get_subhogs_at_level`

        :param hog_id: the query hog id
        :type hog_id: (bytes, str)
        :param str level: the taxonomic level of interest, defaults to None
        :param field: the attribute of the HogLevel table to be returned. Defaults
                      to all attributes of the table.
        """

        if isinstance(hog_id, str):
            hog_id = hog_id.encode("ascii")
        condvars = {"query_fam": self.parse_hog_id(hog_id), "hogid": hog_id}
        query = "(Fam == query_fam) & (ID == hogid)"
        if tab is None:
            tab = self.db.get_node("/HogLevel")
            if level is None:
                query += " & (IsRoot == is_root)"
                condvars["is_root"] = True
            else:
                query += " & (Level == lev)"
                condvars["lev"] = level.encode("ascii")
        try:
            row = next(tab.where(query, condvars=condvars))
            if field is not None:
                if field == "_NROW":
                    return row.nrow
                return row[field]
            return row.fetch_all_fields()
        except StopIteration:
            raise ValueError(
                'HOG-ID/Level combination "{}/{:s}" unknown'.format(
                    hog_id.decode(), level if level is not None else "Root"
                )
            )

    def count_hogs_at_level(self, level):
        """returns the number of HOGs at the requested taxonomic level"""
        try:
            anc_node = self._ancestral_node(level)
            return len(anc_node.Hogs)
        except DBConsistencyError:
            try:
                taxid = self.taxid_from_level(level)
                return len(self.db.get_node("/Hogs_per_Level/tax{}".format(taxid)))
            except tables.NoSuchNodeError:
                return len(self.get_all_hogs_at_level(level))

    def get_all_hogs_at_level(self, level, compare_with=None):
        """returns a :class:`numpy.array` instance with all hogs at the requested level"""
        try:
            anc_node = self._ancestral_node(level)
            hog_data = anc_node.Hogs.read()
        except DBConsistencyError:
            try:
                taxid = self.taxid_from_level(level)
                hog_data = self.db.get_node(
                    "/Hogs_per_Level/tax{}".format(taxid)
                ).read()
            except tables.NoSuchNodeError:
                logger.warning(
                    "Cannot load AncestralGenome nor Hogs_per_Level. extracting from main table. SLOW!"
                )
                hog_data = self.db.get_node("/HogLevel").read_where("Level == level")
        if compare_with is None:
            return hog_data
        compare_hog = self.get_all_hogs_at_level(compare_with)
        diff = compare_levels(compare_hog, hog_data)
        return diff

    def get_roothog_keywords(self, fam) -> str:
        """Retrieve the keywords for a given toplevel HOG.

        :param int fam: numeric id of requested toplevel HOG
        :returns str: keyword of requested toplevel HOG
        """
        try:
            metadata = self.db.get_node("/RootHOG/MetaData")
            buffer = self.db.get_node("/RootHOG/KeywordBuffer")
        except tables.NoSuchNodeError:
            logger.warning("Database has no RootHOG keywords data stored")
            keyword = ""
        else:
            if fam <= 0 or fam > len(metadata):
                raise KeyError("family {} does not exist".format(fam))
            fam_data = metadata[fam - 1]
            if fam_data["FamNr"] != fam:
                logger.warning("/RootHOG/MetaData is not ordered")
                fam_data = metadata.read_where("FamNr == {:d}".format(fam))
            keyword = (
                buffer[
                    fam_data["KeywordOffset"] : fam_data["KeywordOffset"]
                    + fam_data["KeywordLength"]
                ]
                .tobytes()
                .decode()
            )
        return keyword

    def count_hog_members(self, hog_id, level=None):
        """Count the number of members in a (sub)hog.

        If the level is not specified, the deepest level having the given
        hog-id is used.

        :param bytes hog_id: the query hog id
        :param str level: the taxonomic level of interest
        """
        return self.get_hog(hog_id, level, field="NrMemberGenes")

    def get_orthoxml(self, fam, augmented=False):
        """returns the orthoxml of a given toplevel HOG family

        :param fam: numeric id of requested toplevel hog
        :param augmented: boolean flag to indicated whether to return
                          the augmented orthoxml or not. (defaults to not)"""
        idx = read_table_where(self.db.root.OrthoXML.Index, "Fam == {:d}".format(fam))
        if len(idx) < 1:
            raise ValueError("cannot retrieve orthoxml for {}".format(fam))
        idx = idx[0]
        if augmented:
            buf = self.db.root.OrthoXML.BufferAugmented
            rng = slice(
                idx["HogAugmentedBufferOffset"],
                idx["HogAugmentedBufferOffset"] + idx["HogAugmentedBufferLength"],
            )
        else:
            buf = self.db.root.OrthoXML.Buffer
            rng = slice(
                idx["HogBufferOffset"], idx["HogBufferOffset"] + idx["HogBufferLength"]
            )
        return buf[rng].tobytes()

    def get_cached_family_json(self, fam) -> str:
        """returns a json encoded list for protein infos of a family

        This is data contains the protein attributes needed for iHam
        visualization. For older databases that function might raise
        a DBOutdatedError.

        :param fam: numeric id of requested toplevel HOG
        :returns str: json encoded string
        """
        try:
            metadata = self.db.get_node("/RootHOG/MetaData")
            json_buffer = self.db.get_node("/RootHOG/JsonBuffer")
        except tables.NoSuchNodeError as e:
            raise DBOutdatedError(str(e))
        if fam < 0 or fam > len(metadata):
            raise KeyError("family {} does not exist".format(fam))
        fam_data = metadata[fam - 1]
        if fam_data["FamNr"] != fam:
            logger.warning("/RootHOG/MetaData is not ordered")
            fam_data = metadata.read_where("FamNr == {:d}".format(fam))
        json_str = (
            json_buffer[
                fam_data["FamDataJsonOffset"] : fam_data["FamDataJsonOffset"]
                + fam_data["FamDataJsonLength"]
            ]
            .tobytes()
            .decode()
        )
        return json_str

    def _hog_lex_range(self, hog):
        """return the lexographic range of a hog.

        This can be used to search of sub-hogs which are nested in
        the query hog. The semantics is such that
        _hog_lex_range[0] <= hog < _hog_lex_range[1].
        This is equivalent to say that a sub-hog starts with the
        query hog."""
        hog_str = hog.decode() if isinstance(hog, bytes) else hog
        return (
            hog_str.encode("ascii"),
            (hog_str[0:-1] + chr(1 + ord(hog_str[-1]))).encode("ascii"),
        )

    def _ancestral_node(self, level):
        taxid_of_level = self.taxid_from_level(level)
        try:
            return self.db.get_node("/AncestralGenomes/tax{}".format(taxid_of_level))
        except tables.NoSuchNodeError as e:
            logger.exception(
                "cannot find ancestral genome node for taxid:{} ({})".format(
                    taxid_of_level, level
                )
            )
            raise DBConsistencyError(
                f"Ancestral genome node not found: {taxid_of_level}"
            )

    def get_syntentic_hogs(self, level, hog_id, steps=2):
        import warnings

        warnings.warn(
            f"{__name__} is deprecated. please use `get_syntenic_hogs` instead",
            category=DeprecationWarning,
            stacklevel=2,
        )
        return self.get_syntenic_hogs(level, hog_id, steps=steps)

    def get_syntenic_hogs(self, level, hog_id=None, evidence=None, steps=None):
        """Returns a graph of the ancestral synteny

        This method returns a networkx.Graph object with HOGs as nodes
        and weighted edges representing ancestral synteny.

        The evidence is a filter which can be one of `linearized`, `parsimonious` or `any`.
        It defaults to `parsimonious`.

        :param str hog_id: the hog_id of the query HOG
        :param str level: the taxonomic level of the ancestral genome
        :param str evidence: include only edges with evidences up to this level
        :param int steps: number of breadth-first steps to take to get the local
                         neighborhood of the query HOG.
        """
        if self.db_schema_version < (3, 5):
            logger.warning("ancestral synteny: ignoring evidence as outdated database")
            return self._get_syntenic_hogs_obsolete_data(hog_id, level, steps=steps)
        ancestral_node = self._ancestral_node(level)
        if evidence is None:
            evidence = "parsimonious"
        try:
            evidence_enum = ancestral_node.Synteny.get_enum("Evidence")
            evidence = evidence_enum[evidence]
        except KeyError:
            raise ValueError(f"Invalid evidence value {evidence}")
        edge_data = read_table_where(
            ancestral_node.Synteny, "Evidence <= {}".format(evidence)
        )
        edges = (
            (e[0], e[1], {"weight": int(e[2]), "evidence": evidence_enum(e[3])})
            for e in edge_data
        )
        if hog_id is not None:
            hog_row = self.get_hog(hog_id, tab=ancestral_node.Hogs, field="_NROW")
            hogs = ancestral_node.Hogs
        else:
            hog_row = None
            hogs = ancestral_node.Hogs.read()
        return _get_limited_synteny_graph(edges, hogs, hog_row, steps)

    @outdated_database_warning()
    def _get_syntenic_hogs_obsolete_data(self, hog_id, level, steps=2):
        taxid_of_level = self.taxid_from_level(level)
        try:
            edge_data = self.db.get_node(
                "/AncestralSynteny/tax{}".format(taxid_of_level)
            ).read()
        except tables.NoSuchNodeError as e:
            logger.exception(
                "cannot find ancestral synteny data for taxid:{} ({})".format(
                    taxid_of_level, level
                )
            )
            edge_data = []
        edges = ((e[0], e[1], {"weight": e[2], "evidence": "any"}) for e in edge_data)

        hl_tab = self.db.get_node("/HogLevel")
        hogs = hl_tab
        if hog_id is None:
            all_hogs = {}
            for row in hl_tab.where("Level == lev", {"lev": level.encode("utf-8")}):
                all_hogs[row.nrow] = row.fetch_all_fields()
            hog_row = None
            hogs = all_hogs
        else:
            hog_row = self.get_hog(hog_id, level, "_NROW")
        return _get_limited_synteny_graph(edges, hogs, hog_row, steps)

    def taxid_from_level(self, level):
        try:
            taxnodes = self.tax.get_taxnode_from_name_or_taxid(level)
            taxid_of_level = int(taxnodes[0]["NCBITaxonId"])
        except KeyError:
            if level == b"LUCA" or level == "LUCA":
                taxid_of_level = 0
            else:
                raise
        return taxid_of_level

    def get_neighbor_hogs(self, hog_id, level):
        """returns the direct neighbor hogs of a given ancestral hog

        This method returns just the adjacent HOGs from the ancestral
        synteny graph of that level

        :param hog_id: the reference HOG.
        :param level: the taxonomic level of the ancestral genome
        :returns dict: hog_id -> edge_weight of the neighboring hogs
        """
        G = self.get_syntentic_hogs(hog_id=hog_id, level=level, steps=1)
        try:
            neighbors = {nbr: edg["weight"] for nbr, edg in G[hog_id]}
        except KeyError:
            neighbors = {}
        return neighbors

    def get_ancestral_gene_ontology_annotations(
        self,
        level: Union[str, bytes, int],
        hog_id: Union[str, bytes, None] = None,
        as_dataframe: bool = False,
        as_gaf: bool = False,
    ):
        """
        Retrieve the ancestral gene ontology annotations for a hog / all hogs at a taxonomic level

        The method returns the gene ontology annotations stored in the database
        for a given hog_id at a given taxonomic level or for all the hogs at the given
        taxonomic level.

        By default, the result are returned as numpy arrays of type
        :class:`tablefmt.AncestralGeneOntologyTable`, augemtned with a HogID column.
        If `as_dataframe` is set to true, the result will be a pandas dataframe,
        and if `as_gaf` is set to true, a gaf formatted text file with the
        annotations is returned.

        :param level: the taxonomic level for the ancestral data
        :type level: str | bytes | int
        :param hog_id: the hog_id for which the annotations should be retrieved.
                       if not specified, information for all hogs will be used.
        :type hog_id: str | bytes
        :param as_dataframe: return as a pandas dataframe.
        :type as_dataframe: bool
        :param as_gaf: return annotations in GAF format
        :type as_gaf: bool
        """
        ancestral_node = self._ancestral_node(level)
        if hog_id is not None:
            hog = self.get_hog(hog_id, tab=ancestral_node.Hogs)
            query = f"(HogRow == {hog['IdxPerLevelTable']})"
            annots = read_table_where(ancestral_node.GeneOntology, query)
            hog_ids = {hog["IdxPerLevelTable"]: hog["ID"]}
        else:
            annots = ancestral_node.GeneOntology.read()
            hog_ids = ancestral_node.Hogs.read(field="ID")
        annots = numpy.lib.recfunctions.append_fields(
            annots,
            names="HogID",
            data=[hog_ids[row] for row in annots["HogRow"]],
            usemask=False,
        )
        if not as_dataframe and not as_gaf:
            return annots
        # early return if no annotations available
        if len(annots) == 0:
            return "!gaf-version: {}\n".format(GAF_VERSION) if as_gaf else None

        df = pd.DataFrame(annots)
        df["DB"] = "OMA"
        df["DB_Object_ID"] = df["HogID"]
        # 3R DB Object Symbol
        df["DB_Object_Symbol"] = df["DB_Object_ID"]
        # 4O Qualifier
        df["Qualifier"] = ""
        # 5R GO ID
        df["GO_ID"] = df["TermNr"].apply(lambda t: "GO:{:07d}".format(t))
        # 6R DB:Reference
        df["DB:Reference"] = "OMA HogProp"
        # 7R Evidence code
        df["Evidence"] = "IEA"
        # 8O With (or) From
        df["With"] = ""
        # 9R Aspect
        df["Aspect"] = df["GO_ID"].apply(
            lambda t: GOAspect.to_char(self.gene_ontology.term_by_id(t).aspect)
        )
        # 10O DB Object Name
        df["DB_Object_Name"] = ""
        # 11O DB Object Synonym (|Synonym)
        df["Synonym"] = ""
        # 12R DB Object Type
        df["DB_Object_Type"] = "ancestral protein"
        # 13R Taxon (|taxon)
        df["Taxon_ID"] = level
        # 14R Date
        df["Date"] = self.get_conversion_date().strftime("%Y%m%d")
        # 15R Assigned by
        df["Assigned_By"] = df["DB"]
        # 16O Annotation Extension
        df["Annotation_Extension"] = ""
        # 17O Gene Product Form ID
        df["Gene_Product_Form_ID"] = ""

        df = df[GOA.GAF20FIELDS]
        return (
            df
            if not as_gaf
            else (
                "!gaf-version: {}\n".format(GAF_VERSION)
                + "\n".join(df.apply(lambda e: "\t".join(map(str, e)), axis=1))
                + "\n"
            )
        )

    def oma_group_members(self, group_id):
        """get the member entries of an oma group.

        This method returns a numpy array of protein entries that form
        an oma group. If the group id is invalid (not positive
        integer value or a valid Fingerprint), an `InvalidId` Exception
        is raised.

        :param group_id: numeric oma group id or Fingerprint"""
        group_nr = self.resolve_oma_group(group_id)
        pe_tab = self.db.get_node("/Protein/Entries")
        return read_table_where(pe_tab, "OmaGroup == {:d}".format(group_nr))

    def resolve_oma_group(self, group_id):
        if isinstance(group_id, int) and 0 < group_id <= self._nr_oma_groups:
            return group_id
        elif isinstance(group_id, numpy.integer):
            return self.resolve_oma_group(int(group_id))
        elif isinstance(group_id, (bytes, str)):
            if group_id.isdigit():
                return self.resolve_oma_group(int(group_id))
            if isinstance(group_id, str):
                group_id = group_id.encode("utf-8")
            if group_id == b"n/a":
                raise InvalidId("Invalid ID (n/a) for an OMA Group")
            if not self.seq_search.contains_only_valid_chars(group_id):
                raise InvalidId(
                    "Invalid ID: non-amino-accids characters in Fingerprint or sequence pattern"
                )
            if len(group_id) == 7:
                # most likely a fingerprint. let's check that first
                group_meta_tab = self.db.get_node("/OmaGroups/MetaData")
                try:
                    e = next(
                        group_meta_tab.where("(Fingerprint == og)", {"og": group_id})
                    )
                    return int(e["GroupNr"])
                except StopIteration:
                    pass
            # search in suffix array
            entry_nrs = self.seq_search.exact_search(
                group_id.decode(), only_full_length=False
            )
            if len(entry_nrs) == 0:
                raise InvalidId("No sequence contains search pattern")
            group_nrs = {self.entry_by_entry_nr(nr)["OmaGroup"] for nr in entry_nrs}
            group_nrs.discard(0)
            if len(group_nrs) == 1:
                return int(group_nrs.pop())
            elif len(group_nrs) == 0:
                raise InvalidId(
                    "Sequence with pattern '{}' does not belong to any group".format(
                        group_id.decode()
                    )
                )
            else:
                raise AmbiguousID(
                    "sequence pattern matches several oma groups", candidates=group_nrs
                )
        raise InvalidId(
            "Invalid type to determine OMA Group: {} (type: {})".format(
                group_id, type(group_id)
            )
        )

    def oma_group_metadata(self, group_nr):
        """get the metadata associated with an OMA Group

        The metadata contains the fingerprint and the keywords inferred for this group.
        The method returns this information as a dictionary. The parameter must be
        the numeric oma group nr.

        :param int group_nr: a numeric oma group id."""
        if not isinstance(group_nr, (int, numpy.integer)) or group_nr < 0:
            raise InvalidId(
                "Invalid group nr: {} (type: {})".format(group_nr, type(group_nr))
            )
        meta_tab = self.db.get_node("/OmaGroups/MetaData")
        try:
            e = next(meta_tab.where("GroupNr == {:d}".format(group_nr)))
            kw_buf = self.db.get_node("/OmaGroups/KeywordBuffer")
            try:
                grp_size = int(e["NrMembers"])
            except KeyError:
                # fallback if not stored in DB yet
                grp_size = -1

            res = {
                "fingerprint": e["Fingerprint"].decode(),
                "group_nr": int(e["GroupNr"]),
                "keywords": kw_buf[
                    e["KeywordOffset"] : e["KeywordOffset"] + e["KeywordLength"]
                ]
                .tobytes()
                .decode(),
                "size": grp_size,
            }
            return res
        except StopIteration:
            raise InvalidId("invalid group nr")

    def get_nr_oma_groups(self):
        """returns the number of OMA Groups in the database"""
        return self._nr_oma_groups

    @LazyProperty
    def _nr_oma_groups(self):
        tab = self.db.get_node("/Protein/Entries")
        try:
            idx = tab.colindexes["OmaGroup"][-1]
            return int(tab[idx]["OmaGroup"])
        except KeyError:
            hist = self.group_size_histogram("oma")
            return int(hist["Count"].sum())

    def get_nr_toplevel_hogs(self):
        """returns the number of toplevel hogs, i.e. roothogs"""
        hist = self.group_size_histogram("hog")
        return int(hist["Count"].sum())

    def group_size_histogram(self, typ=None):
        """returns a table with two columns, e.g. Size and Count.

        if typ is set to 'oma' or not set, then the data for the
        oma groups is returned. if it is set to 'hog', the data for
        the rootlevel hogs is returned.

        :param typ: either 'oma' or 'hog', defaults to 'oma'"""
        if typ is None or typ.lower() == "oma":
            tabname = "OmaGroup"
        elif typ.lower() == "hog":
            tabname = "OmaHOG"
        else:
            raise ValueError("{} is not a valid group typ".format(typ))
        tab = self.db.get_node("/Summary/{}_size_hist".format(tabname))
        return tab.read()

    def per_species_metadata_retriever(self, genome):
        """return a PerGenomeMetaData instance for a certain query genome to extract
        information on close species based on nr of shared OMA Groups or HOGs"""
        org = self.id_mapper["OMA"].identify_genome(genome)
        return PerGenomeMetaData(self.get_hdf5_handle(), org["UniProtSpeciesCode"])

    def get_sequence(self, entry):
        """get the protein sequence of a given entry as a string

        :param entry: the entry or entry_nr for which the sequence is requested"""
        entry = self.ensure_entry(entry)
        seqArr = self.db.get_node("/Protein/SequenceBuffer")
        seq = seqArr[
            entry["SeqBufferOffset"] : entry["SeqBufferOffset"]
            + entry["SeqBufferLength"]
            - 1
        ]
        return seq.tobytes()

    def get_cdna(self, entry):
        """get the protein sequence of a given entry as a string"""
        entry = self.ensure_entry(entry)
        seqArr = self.db.get_node("/Protein/CDNABuffer")
        seq = seqArr[
            entry["CDNABufferOffset"] : entry["CDNABufferOffset"]
            + entry["CDNABufferLength"]
            - 1
        ]
        return seq.tobytes()

    def get_description(self, entry):
        entry = self.ensure_entry(entry)
        descArr = self.db.get_node("/Protein/DescriptionBuffer")
        desc = descArr[
            entry["DescriptionOffset"] : entry["DescriptionOffset"]
            + entry["DescriptionLength"]
        ]
        return desc.tobytes()

    def get_ec_annotations(self, entry_nr):
        ec_tab = self.db.get_node("/Annotations/EC")
        return [
            row["ECacc"].decode()
            for row in ec_tab.where("EntryNr == {}".format(entry_nr))
        ]

    def get_release_name(self):
        return str(self.db.get_node_attr("/", "oma_version"))

    @LazyProperty
    def release_char(self):
        try:
            release = self.db.get_node_attr("/", "oma_release_char")
        except AttributeError:
            release = ""
        return release

    def get_exons(self, entry_nr):
        genome = (
            self.id_mapper["OMA"]
            .genome_of_entry_nr(entry_nr)["UniProtSpeciesCode"]
            .decode()
        )
        locus_tab = self.db.get_node("/Protein/Locus/{}".format(genome))
        return read_table_where(locus_tab, "EntryNr == {}".format(entry_nr))

    def get_domains(self, entry_nr):
        try:
            return read_table_where(
                self.db.root.Annotations.Domains, "EntryNr == {:d}".format(entry_nr)
            )
        except ValueError as e:
            raise InvalidId("require a numeric entry id, got {}".format(entry_nr))

    def get_representative_entry_of_hog(self, fam):
        """Get the information of the representative entry for a given family (roothog).

        For each family we select a represenative entry that has the most prevalent
        domain architecture. This method returns the entry_nr that we selected, together
        with the domain architecture and its prevalence. In case no representative entry
        has been found, the method raises an :class:`NoReprEntry` Exception.

        :param int fam: The numeric family number."""
        domprev_tab = self.db.get_node("/HOGAnnotations/DomainArchPrevalence")
        try:
            row = next(domprev_tab.where("Fam == {:d}".format(fam)))
            fields = (to_snail_case(z) for z in domprev_tab.dtype.names)
            res = dict(zip(fields, row.fetch_all_fields()))
            res["domains"] = self.get_domains(int(row["ReprEntryNr"]))
            res["prevalence"] = 100.0 * res["prev_count"] / res["fam_size"]
            return res
        except StopIteration:
            raise NoReprEntry()

    def get_prevalent_domains(self, fam):
        # Gets the prevalent domains for a particular top level HOG / family.
        # returns: (family_row, similar_families)
        # family_row contains: family ID, representative entry, DA prevalence.
        # similar_families contains: same, with similarity score. Ordered.
        domprev_tab = self.db.get_node("/HOGAnnotations/DomainArchPrevalence")
        dom2hog_tab = self.db.get_node("/HOGAnnotations/Domains")

        try:
            fam_row = self.get_representative_entry_of_hog(fam)
        except NoReprEntry:
            return None, None

        # Get the family's consensus DA and count them...
        fam_da = collections.Counter(fam_row["domains"]["DomainId"])

        # Retrieve the relevant other families...
        sim_fams = collections.defaultdict(collections.Counter)
        for d in fam_da:
            for hog_with_domain in dom2hog_tab.where("DomainId == {}".format(d)):
                sim_fams[hog_with_domain["Offset"]][d] += 1

        if len(sim_fams) == 0:
            return fam_row, None

        # Now get similar families and order them by similarity
        sim_fams_df = pd.DataFrame(domprev_tab[list(sim_fams.keys())])
        sim_fams_df["sim"] = list(
            map(lambda i: sum((sim_fams[i] & fam_da).values()), sim_fams.keys())
        )

        # Sort by similarity & family size
        sim_fams_df.sort_values(["sim", "FamSize"], inplace=True, ascending=False)
        sim_fams_df.reset_index(drop=True, inplace=True)

        # Prevalence
        sim_fams_df["Prev"] = 100.0 * (
            sim_fams_df["PrevCount"] / sim_fams_df["FamSize"]
        )

        return fam_row, sim_fams_df

    def get_families_with_similar_hog_profile(self, fam, max_nr_similar_fams=50):
        """Retrieves the family nr of families that have a similar
        presence/loss/gain pattern of genes, i.e. potentially
        co-evolving families.

        :param (int,bytes,str) fam: the family of interest, either as a
            hog-id or as a integer identifier

        :param int max_nr_similar_fams: the maximum number of families that
            is returned. Can also be fewer (or even zero).

        :returns dict: a dictionary fam_nr -> np.array indicating which species
            contain at least one gene."""
        if not isinstance(fam, int):
            fam = self.parse_hog_id(fam)
        try:
            if self.hog_profiler is None:
                # create and cache Profiler instance. Not done in constructor
                # as otherwise building of profiles fails.
                self.hog_profiler = Profiler(self)
            return self.hog_profiler.query(fam, k=max_nr_similar_fams)
        except KeyError:
            return None

    def entrynrs_with_ec_annotation(
        self,
        ec: AnyStr,
        limit: Optional[int] = None,
        entrynr_filter: Union[Tuple[int, int], Set[int], None] = None,
    ) -> Set[int]:
        """Search for protein entries with a certain EC annotation.

        The option arguments limit and entrynr_filter can be used to
        limit the amount of entries that are returned and to speed
        up the search time.

        :param ec: EC identifier to be searched
        :type ec: (str, bytes)
        :param limit: maximum number of entries to be returned.
        :type limit: int
        :param entrynr_filter: range (low, high) or set of entry numbers that should be considered
        :type entrynr_filter: Union[Tuple[int, int], Set[int], None]
        :return: Entry numbers of proteins that contain the domain id.
        :rtype: set[int]
        """
        if isinstance(ec, str):
            ec = ec.encode("utf-8")
        query = "(ECacc == ec)"
        ectab = self.get_hdf5_handle().get_node("/Annotations/EC")
        return self._fetch_entry_nrs(
            self._iter_rows_with_entrynr_filter(
                ectab, query, condvars={"ec", ec}, entrynr_filter=entrynr_filter
            ),
            limit,
        )

    def count_ec_annotations(self, ec: AnyStr) -> int:
        """
        Count the number of EC annoations of a certain term in the database.

        :param ec: EC identifier, e.g. "GO:0003676", 50, b"GO:0003674"
        :type ec: (str, bytes)
        :return number of annotations
        :rtype int
        """
        if isinstance(ec, str):
            ec = ec.encode("utf-8")
        ecTab = self.get_hdf5_handle().get_node("/Annotations/EC")
        return count_rows_of_index_column_with_value(ecTab, "ECacc", ec)

    def entrynrs_with_domain_id(
        self,
        domain_id: AnyStr,
        limit: Optional[int] = None,
        entrynr_filter: Union[Tuple[int, int], Set[int], None] = None,
    ) -> Set[int]:
        """Search for protein entries with a certain domain_id.

        The option arguments limit and entrynr_filter can be used to
        limit the amount of entries that are returned and to speed
        up the search time.

        :param domain_id: domain identifier to be searched
        :type domain_id: (str, bytes)
        :param limit: maximum number of entries to be returned.
        :type limit: int
        :param entrynr_filter: range (low, high) or set of entry numbers that should be considered
        :type entrynr_filter:  Union[Tuple[int, int], Set[int], None]
        :return: Entry numbers of proteins that contain the domain id.
        :rtype: set[int]"""

        vars = {
            "dom": domain_id.encode("utf-8")
            if isinstance(domain_id, str)
            else domain_id
        }
        query = "(DomainId == dom)"
        domtab = self.get_hdf5_handle().get_node("/Annotations/Domains")
        return self._fetch_entry_nrs(
            self._iter_rows_with_entrynr_filter(domtab, query, vars, entrynr_filter),
            limit,
        )

    def count_domain_id_annotations(self, domain_id: AnyStr) -> int:
        """
        Count the number of Domain annoations of a certain id in the database.

        :param domain_id: Domain identifier, e.g. "2.60.10.60"
        :type domain_id: (str, bytes)
        :return number of annotations
        :rtype int
        """
        if isinstance(domain_id, str):
            domain_id = domain_id.encode("utf-8")
        ecTab = self.get_hdf5_handle().get_node("/Annotations/Domains")
        return count_rows_of_index_column_with_value(ecTab, "DomainId", domain_id)

    def entrynrs_with_go_annotation(
        self,
        term: Union[AnyStr, int],
        evidence: Optional[AnyStr] = None,
        limit: Optional[int] = None,
        entrynr_filter: Union[Tuple[int, int], Set[int], None] = None,
    ) -> Set[int]:
        """Retrieve protein entry numbers that have a certain GO annotation term

        :param term: numeric term or GO-identifier
        :type term: (int, str, bytes)
        :param evidence: evidence code which is required
        :type evidence: (str, bytes)
        :param limit: maximum number of entries to be returned.
        :type limit: int
        :param entrynr_filter: range (low, high) or set of entry numbers that should be considered
        :type entrynr_filter: Union[Tuple[int, int], Set[int], None]
        :return: Entry numbers of proteins that contain the domain id.
        :rtype: set[int]
        """
        if (isinstance(term, str) and term.startswith("GO:")) or (
            isinstance(term, bytes) and term.startswith(b"GO:")
        ):
            term = term[3:]

        try:
            vars = {"term_nr": int(term)}
        except ValueError:
            raise InvalidId("Invalid GO ID: {}".format(term))

        gotab = self.get_hdf5_handle().get_node("/Annotations/GeneOntology")
        query = "(TermNr == term_nr)"
        if evidence is not None:
            vars["ev"] = evidence.encode("utf-8")
            query += " & (Evidence == ev)"

        return self._fetch_entry_nrs(
            self._iter_rows_with_entrynr_filter(
                gotab, query, condvars=vars, entrynr_filter=entrynr_filter
            ),
            limit,
        )

    def count_go_annotations(self, term: Union[AnyStr, int]) -> int:
        """
        Count the number of GeneOntology annoations of a certain term in the database.
        Annotations are counted multiple time if several evidences or references are
        recorded.

        :param term: numeric term or GO-identifier, e.g. "GO:0003676", 50, b"GO:0003674"
        :type term: (int, str, bytes)
        :return number of annotations
        :rtype int
        """
        if (isinstance(term, str) and term.startswith("GO:")) or (
            isinstance(term, bytes) and term.startswith(b"GO:")
        ):
            term = term[3:]

        try:
            term = int(term)
        except ValueError:
            raise InvalidId("Invalid GO ID: {}".format(term))
        gotab = self.get_hdf5_handle().get_node("/Annotations/GeneOntology")
        return count_rows_of_index_column_with_value(gotab, "TermNr", term)

    def _fetch_entry_nrs(self, row_iter, limit: Optional[int] = None) -> Set[int]:
        entrynrs = set([])
        for row in row_iter:
            entrynrs.add(int(row["EntryNr"]))
            if limit is not None and len(entrynrs) >= limit:
                break
        return entrynrs

    def _iter_rows_with_entrynr_filter(
        self,
        tab: tables.Table,
        query_stub: str,
        condvars: Optional[Mapping] = None,
        entrynr_filter: Union[Tuple[int, int], Set[int], None] = None,
    ) -> Generator:
        query = query_stub
        cond = lambda row: True
        if entrynr_filter is not None:
            rng = min(entrynr_filter), max(entrynr_filter)
            query += " & (EntryNr >= {}) & (EntryNr <= {})".format(*rng)
            if isinstance(entrynr_filter, set):
                cond = lambda row: row["EntryNr"] in entrynr_filter
        for row in tab.where(query, condvars=condvars):
            if cond(row):
                yield row

    def get_gene_ontology_annotations(
        self, entry_nr, stop=None, as_dataframe=False, as_gaf=False
    ):
        """Retrieve the gene ontology annotations for an entry or entry_range

        The method returns the gene ontology annotations stored in the database
        for a given entry (if `stop` parameter is not provided) or for all the
        entries between [entry_nr, stop). Like in slices, the stop entry_nr is
        not inclusive, whereas the entry_nr - the start of the slice - is.

        By default, the result are returned as numpy arrays of type
        :class:`tablefmt.GeneOntologyTable`. If as_dataframe is set to true, the
        result will be a pandas dataframe, and if as_gaf is set to true, a gaf
        formatted text file with the annotations is returned.

        :param int entry_nr: numeric protein entry
        """

        # function to check if an annotation term is obsolete
        def filter_obsolete_terms(term):
            try:
                self.gene_ontology.term_by_id(term)
                return True
            except (KeyError, ValueError):
                return False

        try:
            if stop is None:
                query = "EntryNr == {:d}".format(entry_nr)
            else:
                if not isinstance(stop, (numpy.integer, int)) or stop < entry_nr:
                    raise TypeError(
                        "stop argument needs to be a entry number that is larger than 'entry_nr'"
                    )
                query = "(EntryNr >= {:d}) & (EntryNr < {:d})".format(entry_nr, stop)
            annots = read_table_where(self.db.root.Annotations.GeneOntology, query)

            # for test database we also have some obsolete terms. we need to filter those
            if len(annots) > 0:
                not_obsolete = numpy.vectorize(filter_obsolete_terms)(annots["TermNr"])
                annots = annots[not_obsolete]
        except ValueError as e:
            raise InvalidId("require a numeric entry id, got {}".format(entry_nr))
        if not as_dataframe and not as_gaf:
            return annots

        # early return if no annotations available
        if len(annots) == 0:
            return "!gaf-version: {}\n".format(GAF_VERSION) if as_gaf else None

        df = pd.DataFrame(annots)

        # 1R DB
        df["DB"] = "OMA"
        # 2R DB Object ID
        df["DB_Object_ID"] = df["EntryNr"].apply(self.id_mapper["Oma"].map_entry_nr)
        # 3R DB Object Symbol
        df["DB_Object_Symbol"] = df["DB_Object_ID"]
        # 4O Qualifier
        df["Qualifier"] = ""
        # 5R GO ID
        df["GO_ID"] = df["TermNr"].apply(lambda t: "GO:{:07d}".format(t))
        # 6R DB:Reference
        df["DB:Reference"] = df["Reference"].apply(lambda x: x.decode("ascii"))
        # 7R Evidence code
        df["Evidence"] = df["Evidence"].apply(lambda x: x.decode("ascii"))
        # 8O With (or) From
        df["With"] = ""
        # 9R Aspect
        df["Aspect"] = df["GO_ID"].apply(
            lambda t: GOAspect.to_char(self.gene_ontology.term_by_id(t).aspect)
        )
        # 10O DB Object Name
        df["DB_Object_Name"] = ""
        # 11O DB Object Synonym (|Synonym)
        df["Synonym"] = ""
        # 12R DB Object Type
        df["DB_Object_Type"] = "protein"
        # 13R Taxon (|taxon)
        df["Taxon_ID"] = df["EntryNr"].apply(
            lambda e: "taxon:{:d}".format(
                self.id_mapper["Oma"].genome_of_entry_nr(e)["NCBITaxonId"]
            )
        )
        # 14R Date
        df["Date"] = self.get_conversion_date().strftime("%Y%m%d")
        # 15R Assigned by - TODO: FIX FOR NON OMA!!!
        df["Assigned_By"] = df["DB"]
        # 16O Annotation Extension
        df["Annotation_Extension"] = ""
        # 17O Gene Product Form ID
        df["Gene_Product_Form_ID"] = ""

        df = df[GOA.GAF20FIELDS]
        return (
            df
            if not as_gaf
            else (
                "!gaf-version: {}\n".format(GAF_VERSION)
                + "\n".join(df.apply(lambda e: "\t".join(map(str, e)), axis=1))
                + "\n"
            )
        )

    def get_gene_similarities_hog(self, hog_id):
        """Retrieve the gene similarity matrix for a HOG

        The method returns the 1-D coordinates of the genes of a HOG, indicating
        how close or far they are. Alongside this, the method also returns the
        genes that don't have any go_annotatations present.

        The result is returned as a dictionary

        :param str hog_id: Example 'HOG:0508179'
        """

        hog_members = self.member_of_hog_id(hog_id)

        go_annots_not_fetched = []
        minhash_signatures = {}
        idx_map = {}
        gene_similarity_vals = {}
        total_members = len(hog_members)
        i = 0
        for member in tqdm(
            hog_members,
            disable=len(hog_members) < 500,
            desc="GO-loading {}".format(hog_id),
        ):
            curr_prot_entrynr = member["EntryNr"]
            annos = self.get_gene_ontology_annotations(
                entry_nr=curr_prot_entrynr, as_dataframe=False
            )

            if len(annos) == 0:
                go_annots_not_fetched.append(curr_prot_entrynr)
            else:
                idx_map[i] = curr_prot_entrynr
                h = MinHash()
                for d in annos["TermNr"].astype("unicode"):
                    h.update(d.encode("utf8"))
                minhash_signatures[curr_prot_entrynr] = h
                i += 1

        n = total_members - len(go_annots_not_fetched)
        if n > 0:
            dist_matrix = numpy.zeros((n, n))
            for p1, p2 in tqdm(
                itertools.combinations(range(n), 2),
                disable=n < 50,
                total=n * (n - 1) / 2,
                desc="sim-matrix {}".format(hog_id),
            ):
                score = minhash_signatures[idx_map[p1]].jaccard(
                    minhash_signatures[idx_map[p2]]
                )
                dist_matrix[p1][p2] = 1 - score
                dist_matrix[p2][p1] = 1 - score

            if dist_matrix.max() == 0:
                # all values the same
                positions = numpy.zeros((n, 1))
            else:
                if self.mds is None:
                    from sklearn import manifold

                    self.mds = manifold.MDS(
                        n_components=1,
                        max_iter=100,
                        dissimilarity="precomputed",
                        n_jobs=1,
                    )
                positions = self.mds.fit(dist_matrix).embedding_
            for i in range(len(idx_map)):
                gene_similarity_vals[idx_map[i]] = positions[i][0]
        return go_annots_not_fetched, gene_similarity_vals


class PerGenomeMetaData(object):
    def __init__(self, h5, genome):
        self.h5 = h5
        self.genomes = h5.get_node("/Genome").read(field="UniProtSpeciesCode")
        genome = genome if isinstance(genome, bytes) else genome.encode("ascii")
        try:
            self.genome_idx = numpy.where(self.genomes == genome)[0][0]
        except IndexError:
            raise UnknownSpecies("UniProtSpeciesCode '{}' not known".format(genome))

    def get_most_similar_species(self, limit=None, group_type="OMAGroup", reverse=True):
        overlap_groups = self.h5.get_node(
            "/Summary/shared_{}".format(group_type.lower())
        )[self.genome_idx, :]
        key = numpy.argsort(overlap_groups)
        if reverse:
            limit = limit if limit is None else -limit - 1
            s = slice(None, limit, -1)
        else:
            s = slice(0, limit, None)
        return [
            (self.genomes[k].decode(), int(overlap_groups[k]))
            for k in key[s]
            if k != self.genome_idx
        ]

    def get_least_similar_species(self, **kwargs):
        return self.get_most_similar_species(reverse=False, **kwargs)

    def get_nr_genes_in_group(self, group_type="OMAGroup"):
        return int(
            self.h5.get_node("/Summary/prots_in_{}".format(group_type.lower()))[
                self.genome_idx
            ]
        )


class SequenceSearch(object):
    """
    Contains all the methods for searching the sequence

    TODO: implement taxonomic filtering.
    """

    from .KmerEncoder import DIGITS_AA

    PROTEIN_CHARS = frozenset(map(lambda x: x.decode(), DIGITS_AA))
    PAM100 = pyopa.generate_env(pyopa.load_default_environments()["log_pam1"], 100)

    def __init__(self, db):
        # Backup reference to used DB method.
        self.get_sequence = db.get_sequence

        # Assume the index is stored in the main DB if there is no .idx file
        self.db = db.get_hdf5_handle()
        self.db_idx = (
            self.db
            if not os.path.isfile(self.db.filename + ".idx")
            else tables.open_file(self.db.filename + ".idx", "r")
        )

        # Protein search arrays.
        try:
            self.seq_idx = self.db_idx.root.Protein.SequenceIndex
            if isinstance(self.seq_idx, tables.link.ExternalLink):
                self.seq_idx = self.seq_idx()
            self.kmer_lookup = self.db_idx.root.Protein.KmerLookup
            if isinstance(self.kmer_lookup, tables.link.ExternalLink):
                self.kmer_lookup = self.kmer_lookup()
        except (AttributeError, OSError) as e:
            raise DBConsistencyError(
                "Suffix index for protein sequences is not available: " + str(e)
            )
        self.seq_buff = self.db.root.Protein.SequenceBuffer
        self.n_entries = len(self.db.root.Protein.Entries)

        # Kmer lookup arrays / kmer setup
        self.k = self.kmer_lookup._f_getattr("k")
        self.encoder = KmerEncoder(self.k)
        logger.info("KmerLookup of size k={} loaded".format(self.k))
        self.multienv_align = None

    def get_entry_length(self, ii):
        """Get length of a particular entry."""
        return self.db.root.Protein.Entries[ii - 1]["SeqBufferLength"] - 1

    @LazyProperty
    def entry_idx(self):
        """
        Caches the index lookup part of the SA.
        """
        return self.seq_idx[: self.n_entries]

    def get_entrynr(self, ii):
        """
        Get the entry number(s) corresponding to a location in the sequence
        buffer.
        """
        return numpy.searchsorted(self.entry_idx, ii) + 1

    def contains_only_valid_chars(self, seq):
        """returns true iff `seq` contains only valid AA chars.

        The method ignores the case of the seq, i.e. upper
        or lower case chars both match.

        :param (bytes, str) seq: sequence to be checked
        :returns bool
        """
        if isinstance(seq, bytes):
            seq = seq.decode()
        return all(map(lambda c: c in self.PROTEIN_CHARS, seq.upper()))

    def _sanitise_seq(self, seq):
        """
        Sanitise a string protein sequence. Deletes "invalid" characters.
        TODO: add functionality for biopython sequence / skbio sequence.
        """
        assert type(seq) == str
        return "".join(filter(lambda c: c in self.PROTEIN_CHARS, seq.upper())).encode(
            "ascii"
        )

    def _tax_filter_range(self, en, rng):
        return rng[0] <= en <= rng[1]

    def _tax_filter_set(self, en, taxset):
        return en in taxset

    def search(
        self,
        seq,
        n=None,
        coverage=None,
        is_sanitised=None,
        compute_distance=False,
        entrynr_range=None,
    ):
        """
        Searches the database for entries that match. If can't find an exact
        match performs a kmer + local alignment approach to approximate
        search.
        """
        seq = self._sanitise_seq(seq) if not is_sanitised else seq
        m = self.exact_search(seq, is_sanitised=True, entrynr_range=entrynr_range)
        if len(m) == 0:
            # Do approximate search
            m = self.approx_search(
                seq,
                n=n,
                is_sanitised=True,
                coverage=coverage,
                compute_distance=compute_distance,
                entrynr_range=entrynr_range,
            )
            # TODO: taxonomic filtering.
            return ("approx", m) if m is not [] else None
        else:
            return "exact", m

    def exact_search(
        self, seq, only_full_length=True, is_sanitised=None, entrynr_range=None
    ):
        """
        Performs an exact match search using the suffix array.
        """
        seq = seq if is_sanitised else self._sanitise_seq(seq)
        filt = (
            self._tax_filter_range
            if isinstance(entrynr_range, tuple)
            else self._tax_filter_set
        )
        nn = len(seq)
        if nn > 0:
            z = KeyWrapper(
                self.seq_idx, key=lambda i: self.seq_buff[i : (i + nn)].tobytes()
            )
            ii = bisect_left(z, seq, lo=self.n_entries)

            if ii and ii < len(self.seq_idx) and (z[ii] == seq):
                # Left most found.
                seq_after = seq[:-1] + chr(seq[-1] + 1).encode("utf-8")
                jj = bisect_left(z, seq_after, lo=ii)
                if (jj < len(self.seq_idx) and z[jj] == seq) or z[jj - 1] != seq:
                    raise RuntimeError("suffix index broken. should not happen")

                # Find entry numbers and enr_filter to remove incorrect entries
                return list(
                    filter(
                        lambda e: (
                            ((not only_full_length) or self.get_entry_length(e) == nn)
                            and (entrynr_range is None or filt(e, entrynr_range))
                        ),
                        self.get_entrynr(self.seq_idx[ii:jj]),
                    )
                )

        # Nothing found.
        return []

    def approx_search_no_align(
        self, seq, is_sanitised=False, coverage=0.0, entrynr_range=None
    ):
        """Performs a quick ranking of the best matches based on kmer index.
        The method does not compute an alignment, but just reports the
        entry number with the fraction of kmer matching better than the set coverage

        :param seq: the sequence to be searched
        :type seq: str, bytes
        :param is_sanitised: whether or not the sequence is already sanitised. defaults to false.
        :type is_sanitised: bool
        :param coverage: the minimum fraction of covered kmers by the target sequence
        :type coverage: float
        :param entrynr_range: target entry number range as a tuple (min, max) or set for filtering
        :type entrynr_range: set[int], tuple[int, int]
        :returns: A list of tuples with (entry_nr, fraction_of_matched_kmers)
        """
        seq = seq if is_sanitised else self._sanitise_seq(seq)
        tax_filt = (
            self._tax_filter_range
            if isinstance(entrynr_range, tuple)
            else self._tax_filter_set
        )

        # 1. Do kmer counting vs entry numbers TODO: switch to np.unique?
        c = collections.Counter()
        for z in map(
            lambda kmer: numpy.unique(self.kmer_lookup[int(kmer)]),
            self.encoder.decompose(seq),
        ):
            c.update(z)

        # 2. Filter to top n if necessary
        z = len(seq) - self.k + 1
        cut_off = coverage * z
        entries = [
            (enr, (cnts / z))
            for enr, cnts in c.items()
            if cnts >= cut_off
            and (entrynr_range is None or tax_filt(enr, entrynr_range))
        ]
        return entries

    def approx_search(
        self,
        seq,
        n=None,
        is_sanitised=None,
        coverage=None,
        compute_distance=False,
        entrynr_range=None,
        return_kmer_hits=False,
    ):
        """
        Performs an approximate match search using the kmer index.
        For the n best matching candidates based on kmer lookup, a
        Smith-Waterman alignment is computed.

        If compute_distance is set to True, an ML distance estimate
        of the alignment is returned as part of the resulting dictionary

        :param seq: the query sequence to be searched
        :type seq: str, bytes
        :param int n: number of maximum returned entries that match query
        :param bool is_sanitised: whether or not the sequence is already sanitised. defaults to false.
        :param float coverage: the minimum fraction of covered kmers by the target sequence
        :param entrynr_range: target entry number range as a tuple (min, max) or set for filtering
        :type entrynr_range: set[int], tuple[int, int]
        :param bool return_kmer_hits: whether or not the full list of matched kmer entries should be
        returned. if set to True, the return value will be a tuple instead of a single list

        :returns: list of matched entries, each element is a tuple with the entry_nr and a dictionary
        containing the score, alignment and distance estimates.
        If return_kmer_hits is set to True, also the full list of matching entries is returned.
        """
        seq = seq if is_sanitised else self._sanitise_seq(seq)
        n = n if n is not None else 50
        coverage = 0.0 if coverage is None else coverage

        kmer_hits = self.approx_search_no_align(
            seq, is_sanitised=True, coverage=coverage, entrynr_range=entrynr_range
        )
        c = sorted(kmer_hits, reverse=True, key=lambda x: x[1])
        if n > 0:
            c = c[:n]

        # 3. Do local alignments and return count / score / alignment
        res = []
        if len(c) > 0:
            res = sorted(
                [
                    (
                        m[0],
                        {
                            "kmer_coverage": m[1],
                            "score": a[0],
                            "alignment": a[1],
                            "distance": a[2] if compute_distance else None,
                            "distvar": a[3] if compute_distance else None,
                        },
                    )
                    for (m, a) in self._align_entries(seq, c, compute_distance)
                ],
                key=lambda q: q[1]["score"],
                reverse=True,
            )
        if return_kmer_hits:
            return res, kmer_hits
        else:
            return res

    def _align_entries(self, seq, matches, compute_distance=False):
        # Does the alignment for the approximate search
        def align(s1, s2s, env, aligned):
            for s2 in s2s:
                z = pyopa.align_double(s1, s2, env, False, False, True)
                a = pyopa.align_strings(s1, s2, env, False, z)
                if compute_distance:
                    score, pam, pamvar = self.multienv_align.estimate_pam(*a[0:2])
                    res_ds = (
                        score,
                        (
                            (a[0].convert_readable(), (z[3], z[1])),
                            (a[1].convert_readable(), (z[4], z[2])),
                        ),
                        pam,
                        pamvar,
                    )
                else:
                    res_ds = (
                        z[0],
                        (
                            (a[0].convert_readable(), (z[3], z[1])),
                            (a[1].convert_readable(), (z[4], z[2])),
                        ),
                    )
                aligned.append(res_ds)

        if compute_distance and self.multienv_align is None:
            # lazy loading of MultipleAlEnv datastructure
            envs = pyopa.load_default_environments()
            self.multienv_align = pyopa.MutipleAlEnv(
                envs["environments"], envs["log_pam1"]
            )

        aligned = []
        query = pyopa.Sequence(seq.decode("ascii"))
        entries = list(
            map(
                lambda m: pyopa.Sequence(self.get_sequence(int(m[0])).decode("ascii")),
                matches,
            )
        )
        t = threading.Thread(target=align, args=(query, entries, self.PAM100, aligned))
        t.start()
        t.join()
        assert len(aligned) > 0, "Alignment thread crashed."
        return zip(matches, aligned)


class OmaIdMapper(object):
    def __init__(self, db):
        self.genome_table = db.get_hdf5_handle().root.Genome.read()
        self._entry_off_keys = self.genome_table.argsort(order=("EntryOff"))
        self._genome_keys = self.genome_table.argsort(order=("UniProtSpeciesCode"))
        self._taxid_keys = self.genome_table.argsort(order=("NCBITaxonId"))
        self._sciname_keys = self.genome_table.argsort(order=("SciName"))
        self._omaid_re = re.compile(r"(?P<genome>[A-Z][A-Z0-9]{4})(?P<nr>\d+)")
        self._db = db
        self._approx_genome_matcher = self._init_fuzzy_matcher_with_genome_infos()

    def _init_fuzzy_matcher_with_genome_infos(self):
        values = []
        maps_to = []
        fields = {
            "SciName",
            "SynName",
            "CommonName",
            "UniProtSpeciesCode",
        }.intersection(self.genome_table.dtype.names)
        for col in fields:
            for row in range(len(self.genome_table)):
                val = self.genome_table[col][row].decode()
                if len(val) > 0:
                    values.append(val)
                    maps_to.append(row)
        return FuzzyMatcher(values, maps_to, rel_sim_cutoff=0.6)

    def genome_of_entry_nr(self, e_nr):
        """returns the genome code belonging to a given entry_nr"""
        idx = self.genome_table["EntryOff"].searchsorted(
            e_nr - 1, side="right", sorter=self._entry_off_keys
        )
        return self.genome_table[self._entry_off_keys[idx - 1]]

    def map_entry_nr(self, entry_nr):
        genome = self.genome_of_entry_nr(entry_nr)
        return "{0:s}{1:05d}".format(
            genome["UniProtSpeciesCode"].decode(), entry_nr - genome["EntryOff"]
        )

    def genome_from_UniProtCode(self, code):
        code = code.encode("ascii")
        idx = self.genome_table["UniProtSpeciesCode"].searchsorted(
            code, sorter=self._genome_keys
        )
        try:
            genome = self.genome_table[self._genome_keys[idx]]
        except IndexError:
            raise UnknownSpecies("{} is unknown".format(code))

        if genome["UniProtSpeciesCode"] != code:
            raise UnknownSpecies("{} is unknown".format(code))
        return genome

    def genome_from_SciName(self, name):
        name = name.encode("ascii")
        idx = self.genome_table["SciName"].searchsorted(name, sorter=self._sciname_keys)
        try:
            genome = self.genome_table[self._sciname_keys[idx]]
        except IndexError:
            raise UnknownSpecies("{} is unknown".format(name))

        if genome["SciName"] != name:
            raise UnknownSpecies("{} is unknown".format(name))
        return genome

    def genome_from_taxid(self, taxid):
        try:
            taxid = int(taxid)
            idx = self.genome_table["NCBITaxonId"].searchsorted(
                taxid, sorter=self._taxid_keys
            )
            genome = self.genome_table[self._taxid_keys[idx]]
        except (IndexError, ValueError):
            raise UnknownSpecies('TaxonId "{}" is unknown'.format(taxid))
        if genome["NCBITaxonId"] != taxid:
            raise UnknownSpecies('TaxonId "{}" is unknown'.format(taxid))
        return genome

    def identify_genome(self, code):
        """identify genome based on either a UniProtSpeciesCode,
        NCBI Taxonomy Id or species name"""
        try:
            code = int(code)
        except ValueError:
            if len(code) == 5:
                try:
                    return self.genome_from_UniProtCode(code.upper())
                except UnknownSpecies:
                    pass
            return self.genome_from_SciName(code)
        return self.genome_from_taxid(code)

    def approx_search_genomes(self, pattern, scores=False):
        candidates = self._approx_genome_matcher.search_approx(pattern)
        genomes = [Genome(self._db, self.genome_table[z[2]]) for z in candidates]
        if scores:
            return genomes, [z[0] for z in candidates]
        else:
            return genomes

    def omaid_to_entry_nr(self, omaid):
        """returns the internal numeric entrynr from a
        UniProtSpeciesCode+nr id. this is the inverse
        function of 'map_entry_nr'."""
        match = self._omaid_re.match(omaid)
        if match is None:
            raise InvalidOmaId(omaid)
        code, nr = match.group("genome"), int(match.group("nr"))
        genome = self.genome_from_UniProtCode(code)
        if nr <= 0 or nr > genome["TotEntries"]:
            raise InvalidOmaId(omaid)
        return genome["EntryOff"] + int(match.group("nr"))

    def genome_range(self, query):
        """returns the internal range of EntryNr associated with
        'query'. 'query' can be either a numeric id of a protein
        or a UniProtSpeciesCode of a genome. If 'query' is unknown
        by the database, an InvalidOmaId exception is raised.

        The return range is a tuple of length two, and the numbers
        indicated the *inclusive* boundaries, e.g. (1,5) indicates
        that the entries 1,2,3,4 and 5 belong to the query species"""
        if isinstance(query, (int, numpy.integer)):
            genome_row = self.genome_of_entry_nr(query)
            if query <= 0 or query > genome_row["EntryOff"] + genome_row["TotEntries"]:
                raise InvalidOmaId(query)
        else:
            genome_row = self.genome_from_UniProtCode(query)
        return (
            genome_row["EntryOff"] + 1,
            genome_row["EntryOff"] + genome_row["TotEntries"],
        )

    def species_ordering(self, root=None):
        """get ordering of the genomes with respect to taxonomy.

        This method returns a linear ordering of all the genomes with
        respect to their lineage, i.e. genomes that are evolutionary
        "close" to each other appear close in the ordering.
        Optionally, one can give a root genome, that will be the species
        the ordering is going to start with.

        :param root: UniProtSpeciesCode of the root genome.
        :returns: a list of species codes in the correct order."""
        if root is None:
            root = self.genome_table[0]["UniProtSpeciesCode"]
        root_genome = self.genome_from_UniProtCode(root)
        lins = {
            g["UniProtSpeciesCode"]: [
                lev["Name"] for lev in self._db.tax.get_parent_taxa(g["NCBITaxonId"])
            ][::-1]
            for g in self.genome_table
        }
        root_lin = lins[root_genome["UniProtSpeciesCode"]]
        sort_key = {}
        for g, lin_g in lins.items():
            for k in range(min(len(root_lin), len(lin_g))):
                if root_lin[k] != lin_g[k]:
                    k -= 1
                    break
            sort_key[g] = (-k, lin_g)
        sorted_genomes = sorted(list(sort_key.keys()), key=lambda g: sort_key[g])
        return {g.decode(): v for v, g in enumerate(sorted_genomes)}


class HogIdForwardMapper(object):
    def __init__(self, db: Database, hogmap_filename=None):
        if hogmap_filename is None:
            hogmap_filename = db.get_hdf5_handle().filename.replace(".h5", ".hogmap.h5")
        self.h5_hogmap = tables.open_file(hogmap_filename, mode="r")
        self.db = db
        self.db.register_on_close(self.close)

    def close(self):
        self.h5_hogmap.close()
        self.db.unregister_on_close(self.close)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def map_hogid(self, hogid):
        hogid = hogid.encode("utf-8") if isinstance(hogid, str) else hogid
        cand_iter = self.h5_hogmap.get_node("/hogmap").where(
            "Old == hogid",
            {"hogid": hogid},
        )
        return {c["New"].decode(): float(c["Jaccard"]) for c in cand_iter}


class FuzzyMatcher(object):
    def __init__(
        self,
        values,
        maps_to=None,
        rel_sim_cutoff=0.8,
        gram_size_lower=3,
        gram_size_upper=5,
    ):
        """FuzzyMatcher allows to search for approximate matches of a list of values.
        It is a thin wrapper to the :class:`fuzzyset.FuzzySet datastructure.

        FuzzyMatcher can be initialized with either a pure list of values or including
        a seperate list with mapping objects. The values (repetitions are possible)
        indicate to what object they should be mapped
        On a search (see :meth:`search_approx`) the object associated with the key
        will be returned.

        :param values: an iterable/mapping
        """
        if maps_to is not None:
            self.fuzzySet = fuzzyset.FuzzySet(
                rel_sim_cutoff=rel_sim_cutoff,
                gram_size_lower=gram_size_lower,
                gram_size_upper=gram_size_upper,
            )
            self.mapping = collections.defaultdict(list)
            for val, map_source in zip(values, maps_to):
                self.fuzzySet.add(val)
                self.mapping[val].append(map_source)
        else:
            self.fuzzySet = fuzzyset.FuzzySet(
                values,
                rel_sim_cutoff=rel_sim_cutoff,
                gram_size_lower=gram_size_lower,
                gram_size_upper=gram_size_upper,
            )
            self.mapping = None

    def search_approx(self, key):
        matches = self.fuzzySet.get(key, [])
        if self.mapping:
            bests = {}
            for score, val in matches:
                sources = self.mapping[val]
                for src in sources:
                    if src not in bests or score > bests[src][0]:
                        bests[src] = (score, val, src)
            matches = list(bests.values())
        return matches


class IDResolver(object):
    def __init__(self, db):
        entry_nr_col = db.get_hdf5_handle().root.Protein.Entries.cols.EntryNr
        self.max_entry_nr = entry_nr_col[int(entry_nr_col.index[-1])]
        self._db = db

    def _from_numeric(self, e_id):
        nr = int(e_id)
        if not 0 < nr <= self.max_entry_nr:
            raise InvalidId("{0:d} out of protein range: {1:}".format(nr, e_id))
        return nr

    def _from_omaid(self, e_id):
        return int(self._db.id_mapper["OMA"].omaid_to_entry_nr(e_id))

    def search_xrefs(self, e_id, return_seq_modified=False):
        """search for all xrefs. TODO: what happens if xref is ambiguous?"""
        xref_res = self._db.id_mapper["XRef"].search_xref(e_id)
        res = set([x["EntryNr"] for x in xref_res])
        if len(res) == 0:
            # let's try to mach as substring using suffix array case insensitive
            xref_res = self._db.id_mapper["XRef"].search_xref(
                e_id, match_any_substring=True
            )
            res = set([x["EntryNr"] for x in xref_res])
            if len(res) == 0:
                raise InvalidId(e_id)
        if len(res) > 1:
            # check whether its only different isoforms, then return canonical isoform
            splice_variants = set(
                [
                    x["AltSpliceVariant"]
                    for x in (self._db.entry_by_entry_nr(eNr) for eNr in res)
                ]
            )
            logger.info(
                "xref {} has {} entries, {} splice variants".format(
                    e_id, len(res), len(splice_variants)
                )
            )
            if len(splice_variants) > 1 or 0 in splice_variants:
                raise AmbiguousID('Cross-ref "{}" is ambiguous'.format(e_id), res)
            else:
                res = splice_variants
        return_val = int(res.pop())
        if return_seq_modified:
            idx = numpy.where(xref_res["EntryNr"] == return_val)
            return return_val, bool(4 in xref_res[idx]["Verification"])
        return return_val

    def resolve(self, e_id, check_if_modified=False):
        """maps an id to the entry_nr of the current OMA release."""
        try:
            nr = self._from_numeric(e_id)
        except ValueError:
            try:
                nr = self._from_omaid(e_id)
            except (InvalidOmaId, UnknownSpecies) as e:
                nr = self.search_xrefs(e_id, return_seq_modified=check_if_modified)
        if check_if_modified and not isinstance(nr, tuple):
            nr = (nr, False)
        return nr

    @timethis(logging.DEBUG)
    def search_protein(self, query: str, limit=None, entrynr_range=None):
        candidates = collections.defaultdict(dict)
        try:
            nr = self._from_numeric(query)
            if entrynr_range is None or entrynr_range[0] <= nr <= entrynr_range[1]:
                candidates[nr]["numeric_id"] = [query]
        except ValueError:
            pass
        try:
            nr = self._from_omaid(query)
            if entrynr_range is None or entrynr_range[0] <= nr <= entrynr_range[1]:
                candidates[nr]["omaid"] = [query]
        except (InvalidOmaId, UnknownSpecies) as e:
            pass
        id_res = self._db.id_mapper["XRef"].search_id(
            query, limit=limit, entrynr_range=entrynr_range
        )
        for nr, res_dict in id_res.items():
            candidates[nr].update(res_dict)
        if len(id_res) < 5 and self._db.desc_searcher is not None:
            try:
                desc_res = self._db.desc_searcher.search_term(query, limit)
                desc_res.sort()
                for row_nr in desc_res[:limit]:
                    nr = row_nr + 1
                    if (
                        entrynr_range is None
                        or entrynr_range[0] <= nr <= entrynr_range[1]
                    ):
                        candidates[nr]["Description"] = [
                            self._db.get_description(nr).decode()
                        ]
            except SuffixIndexError as e:
                logger.warning(e)
                logger.warning("No Descriptions searched")
        return candidates


class Taxonomy(object):
    """Taxonomy provides an interface to navigate the taxonomy data.

    The input data is the same as what is stored in the Database in
    table "/Taxonomy"."""

    def __init__(self, data, genomes=None, _valid_levels=None):
        if not isinstance(data, numpy.ndarray):
            raise ValueError("Taxonomy expects a numpy table.")
        self.genomes = genomes if genomes is not None else {}
        self.tax_table = data
        self.taxid_key = self.tax_table.argsort(order=("NCBITaxonId"))
        self.parent_key = self.tax_table.argsort(order=("ParentTaxonId"))
        self._add_luca_if_needed()
        self.all_hog_levels = _valid_levels
        self._approx_matcher = None
        if _valid_levels is None:
            self._load_valid_taxlevels()

    def _load_valid_taxlevels(self):
        forbidden_chars = re.compile(r"[^A-Za-z0-9()/_,.: -]")
        try:
            with open(os.environ["DARWIN_BROWSERDATA_PATH"] + "/TaxLevels.drw") as f:
                taxStr = f.read()
            tax_json = json.loads(("[" + taxStr[14:-3] + "]").replace("'", '"'))
            self.all_hog_levels = frozenset(
                [
                    t.encode("ascii")
                    for t in tax_json
                    if forbidden_chars.search(t) is None
                ]
            )
        except (IOError, KeyError):
            self.all_hog_levels = frozenset(
                [
                    l
                    for l in self.tax_table["Name"]
                    if forbidden_chars.search(l.decode()) is None
                ]
            )

    def _table_idx_from_numeric(self, tids):
        i = self.tax_table["NCBITaxonId"].searchsorted(tids, sorter=self.taxid_key)
        idx = self.taxid_key[i]
        if (self.tax_table[idx]["NCBITaxonId"] != tids).any():
            unkn = tids[self.tax_table[idx]["NCBITaxonId"] != tids]
            raise InvalidTaxonId("{} are invalid/unknown taxonomy ids".format(unkn))
        return idx

    def _get_root_taxon(self):
        i1 = self.tax_table["ParentTaxonId"].searchsorted(0, sorter=self.parent_key)
        i2 = self.tax_table["ParentTaxonId"].searchsorted(
            0, sorter=self.parent_key, side="right"
        )
        if i2 - i1 == 0:
            raise DBConsistencyError(
                "Not a single root in Taxonomy: {}".format(
                    self.tax_table[self.parent_key[i1]]
                )
            )
        elif i2 - i1 == 1:
            res = self.tax_table[self.parent_key[i1]]
        else:
            res = numpy.array([(0, -1, b"LUCA")], dtype=self.tax_table.dtype)[0]
        return res

    def _add_luca_if_needed(self):
        i1 = self.tax_table["ParentTaxonId"].searchsorted(0, sorter=self.parent_key)
        i2 = self.tax_table["ParentTaxonId"].searchsorted(
            0, sorter=self.parent_key, side="right"
        )
        if i2 - i1 > 1:
            self.tax_table = numpy.append(
                self.tax_table,
                numpy.array([(0, -1, b"LUCA")], dtype=self.tax_table.dtype),
            )
            self.taxid_key = self.tax_table.argsort(order=("NCBITaxonId"))
            self.parent_key = self.tax_table.argsort(order=("ParentTaxonId"))

    def _taxon_from_numeric(self, tid):
        idx = self._table_idx_from_numeric(tid)
        return self.tax_table[idx]

    def _direct_children_taxa(self, tid):
        i = self.tax_table["ParentTaxonId"].searchsorted(tid, sorter=self.parent_key)
        idx = []
        while (
            i < len(self.parent_key)
            and self.tax_table[self.parent_key[i]]["ParentTaxonId"] == tid
        ):
            idx.append(self.parent_key[i])
            i += 1
        return self.tax_table.take(idx)

    def approx_search(self, pattern):
        if self._approx_matcher is None:
            valid_levels_as_str = (l.decode() for l in self.all_hog_levels)
            self._approx_matcher = FuzzyMatcher(values=valid_levels_as_str)
        candidates = self._approx_matcher.search_approx(pattern)
        return candidates

    def get_parent_taxa(self, query, root=0):
        """Get array of taxonomy entries leading towards the
        root of the taxonomy.

        :param query: the starting taxonomy level
        :param root: parameter to identify the root-taxid. defaults to 0.
                     if set to -1, LUCA will be included into the list.
        """
        idx = []
        parent = query
        count = 0
        while parent != root:
            i = self._table_idx_from_numeric(parent)
            idx.append(i)
            tmp = self.tax_table[i]["ParentTaxonId"]
            if tmp == parent:
                raise InvalidTaxonId("{0:d} has itself as parent".format(tmp))
            parent = tmp
            count += 1
            if count > 100:
                raise InvalidTaxonId(
                    "{0:d} exceeds max depth of 100. Infinite recursion?".format(query)
                )
        return self.tax_table.take(idx)

    def _get_taxids_from_any(self, it, skip_missing=True):
        if not isinstance(it, numpy.ndarray):
            try:
                it = numpy.fromiter(it, dtype="i4")
            except ValueError:
                it = numpy.fromiter(it, dtype="S255")
        if it.dtype.type is numpy.string_:
            try:
                ns = self.name_key
            except AttributeError:
                ns = self.name_key = self.tax_table.argsort(order="Name")
            idxs = self.tax_table["Name"].searchsorted(it, sorter=ns)
            idxs = numpy.clip(idxs, 0, len(ns) - 1)
            taxs = self.tax_table[ns[idxs]]
            keep = taxs["Name"] == it
            if not skip_missing and not keep.all():
                raise KeyError("not all taxonomy names could be found")
            res = taxs["NCBITaxonId"][keep]
        else:
            res = it
        return res

    def get_subtaxonomy_rooted_at(self, root, collapse=True):
        rid = self._get_taxids_from_any([root])
        subtree = [rid]

        def get_children(id):
            children = self._direct_children_taxa(id)
            if len(children) > 0:
                for child in children:
                    child_id = child["NCBITaxonId"]
                    subtree.append(child_id)
                    get_children(child_id)

        get_children(rid)
        return self.get_induced_taxonomy(subtree, collapse=collapse)

    def get_taxnode_from_name_or_taxid(self, query):
        if isinstance(query, (bytes, str, int)):
            query = [query]
        tids = self._get_taxids_from_any(query, skip_missing=False)
        return self._taxon_from_numeric(tids)

    def get_taxid_of_extent_genomes(self):
        """returns a list of ncbi taxon ids of the extent genomes within the taxonomy"""

        def _traverse(node):
            children = self._direct_children_taxa(node["NCBITaxonId"])
            for child in children:
                _traverse(child)

            if len(children) == 0 or (int(node["NCBITaxonId"]) in self.genomes):
                extent_genomes.append(int(node["NCBITaxonId"]))

        extent_genomes = []
        _traverse(self._get_root_taxon())
        return extent_genomes

    def get_induced_taxonomy(self, members, collapse=True, augment_parents=False):
        """Extract the taxonomy induced by a given set of `members`.

        This method allows to extract the part which is induced by a
        given set of levels and leaves that should be part of the
        new taxonomy. `members` must be an iterable, the levels
        must be either numeric taxids or scientific names.

        Unless `augment_parents` is set to true, the resulting sub-taxonomy
        will only contain levels that are specified in `members`. If
        `augment_parents` is set to True, also all parent nodes of the
        levels passed in members are considered for the sub-taxonomy.

        :param iter members: an iterable containing the levels
            and leaves that should remain in the new taxonomy. can be
            either axonomic ids or scientific names.

        :param bool collapse: whether or not levels with only one child
            should be skipped or not. This defaults to True

        :param bool augment_parents: whether or not to consider parent
            levels of members for the resulting taxonomy."""

        taxids_to_keep = numpy.sort(self._get_taxids_from_any(members))
        if augment_parents:
            # find all the parents of all the members, add them to taxids_to_keep
            additional_levels = set([])
            for cur_tax in taxids_to_keep:
                try:
                    additional_levels.update(
                        set(self.get_parent_taxa(cur_tax)["NCBITaxonId"])
                    )
                except KeyError:
                    logger.info("{} seems not to exist in Taxonomy".format(cur_tax))
                    pass
            # add and remove duplicates
            all_levels = numpy.append(taxids_to_keep, list(additional_levels))
            taxids_to_keep = numpy.unique(all_levels)

        idxs = numpy.searchsorted(
            self.tax_table["NCBITaxonId"], taxids_to_keep, sorter=self.taxid_key
        )
        idxs = numpy.clip(idxs, 0, len(self.taxid_key) - 1)
        subtaxdata = self.tax_table[self.taxid_key[idxs]]
        if not numpy.alltrue(subtaxdata["NCBITaxonId"] == taxids_to_keep):
            raise KeyError("not all levels in members exists in this taxonomy")

        updated_parent = numpy.zeros(len(subtaxdata), "bool")
        for i, cur_tax in enumerate(taxids_to_keep):
            if updated_parent[i]:
                continue
            # get all the parents and check which ones we keep in the new taxonomy.
            parents = self.get_parent_taxa(cur_tax)["NCBITaxonId"]
            mask = numpy.in1d(parents, taxids_to_keep)
            # find the position of them in subtaxdata (note: subtaxdata and
            # taxids_to_keep have the same ordering).
            new_idx = taxids_to_keep.searchsorted(parents[mask])
            taxids = taxids_to_keep[new_idx]
            # parent taxid are ncbitaxonids shifted by one position!
            parents = numpy.roll(taxids, -1)
            parents[-1] = 0
            subtaxdata["ParentTaxonId"][new_idx] = parents
            updated_parent[new_idx] = True

        if collapse:
            nr_children = collections.defaultdict(int)
            for p in subtaxdata["ParentTaxonId"]:
                nr_children[p] += 1
            rem = [
                p
                for (p, cnt) in nr_children.items()
                if cnt == 1 and p != 0 and p not in self.genomes
            ]
            if len(rem) > 0:
                idx = taxids_to_keep.searchsorted(rem)
                return self.get_induced_taxonomy(numpy.delete(taxids_to_keep, idx))
        return Taxonomy(
            subtaxdata, genomes=self.genomes, _valid_levels=self.all_hog_levels
        )

    def newick(self):
        """Get a Newick representation of the Taxonomy

        Note: as many newick parsers do not support quoted labels,
        the method instead replaces spaces with underscores."""

        def newick_enc(s):
            return s.translate({ord(" "): "_", ord("("): "[", ord(")"): "]"})

        def _rec_newick(node):
            children = []
            for child in self._direct_children_taxa(node["NCBITaxonId"]):
                children.append(_rec_newick(child))

            if len(children) == 0:
                return newick_enc(node["Name"].decode())
            else:
                if int(node["NCBITaxonId"]) in self.genomes:
                    # this is a special case where the current internal level is also
                    # an extant species in OMA. we resolve this by adding the current
                    # level also as an extra child
                    children.append(
                        newick_enc(
                            "{:s} (disambiguate {:s})".format(
                                node["Name"].decode(),
                                self.genomes[
                                    int(node["NCBITaxonId"])
                                ].uniprot_species_code,
                            )
                        )
                    )

                t = ",".join(children)
                return "(" + t + ")" + newick_enc(node["Name"].decode())

        return _rec_newick(self._get_root_taxon()) + ";"

    def as_dict(self):
        """Encode the Taxonomy as a nested dict.

        This representation can for example be used to serialize
        a Taxonomy in json format."""

        def _rec_phylogeny(node):
            res = {"name": node["Name"].decode(), "id": int(node["NCBITaxonId"])}
            children = []
            for child in self._direct_children_taxa(node["NCBITaxonId"]):
                children.append(_rec_phylogeny(child))
            if len(children) > 0:
                if res["id"] in self.genomes:
                    # this is a special case where the current internal level is also
                    # an extant species in OMA. we resolve this by adding the current
                    # level also as an extra child
                    code = self.genomes[res["id"]].uniprot_species_code
                    node_cpy = {
                        "name": "{:s} (disambiguate {:s})".format(res["name"], code),
                        "id": res["id"],
                        "code": code,
                    }
                    children.append(node_cpy)
                res["children"] = children
            else:
                try:
                    g = self.genomes[res["id"]]
                    res["code"] = g.uniprot_species_code
                except KeyError:
                    pass
            return res

        return _rec_phylogeny(self._get_root_taxon())

    def as_phyloxml(self):
        """Encode the Taxonomy as phyloxml output"""

        def _rec_phyloxml(node):
            n = et.Element("clade")
            tax = et.SubElement(n, "taxonomy")
            id_ = et.SubElement(tax, "id", provider="uniprot")
            id_.text = str(node["NCBITaxonId"])

            children = []
            for child in self._direct_children_taxa(node["NCBITaxonId"]):
                children.append(_rec_phyloxml(child))
            if len(children) > 0 and int(node["NCBITaxonId"]) in self.genomes:
                # this is a special case where the current internal level is also
                # an extant species in OMA. we resolve this by adding the current
                # level also as an extra child
                cp_n = et.Element("clade")
                cp_tax = et.SubElement(cp_n, "taxonomy")
                cp_id = et.SubElement(cp_tax, "id", provider="uniprot")
                cp_id.text = str(node["NCBITaxonId"])
                cp_code = et.SubElement(cp_tax, "code")
                cp_code.text = self.genomes[
                    int(node["NCBITaxonId"])
                ].uniprot_species_code
                cp_sci = et.SubElement(cp_tax, "scientific_name")
                cp_sci.text = "{:s} (disambiguate {:s})".format(
                    node["Name"].decode(), cp_code.text
                )
                children.append(cp_n)

            if len(children) == 0:
                try:
                    g = self.genomes[int(node["NCBITaxonId"])]
                    code = et.SubElement(tax, "code")
                    code.text = g.uniprot_species_code
                except (ValueError, KeyError):
                    pass
            sci = et.SubElement(tax, "scientific_name")
            sci.text = node["Name"].decode()
            n.extend(children)
            return n

        root = et.Element("phyloxml", xmlns="http://www.phyloxml.org")
        phylo = et.SubElement(root, "phylogeny", rooted="true", rerootable="false")
        name = et.SubElement(phylo, "name")
        name.text = "(Partial) species phylogeny from OMA Browser"
        phylo.append(_rec_phyloxml(self._get_root_taxon()))

        return et.tostring(root, encoding="utf-8")


class IdMapperFactory(object):
    def __init__(self, db_obj):
        self.db = db_obj
        self.mappers = {}

    def __getitem__(self, idtype):
        return self.get_mapper(idtype)

    def get_mapper(self, idtype):
        try:
            mapper = self.mappers[idtype]
        except KeyError:
            mapper = None
            for inst in (str(idtype).title() + "IdMapper", str(idtype)):
                try:
                    mapper = globals()[inst](self.db)
                    self.mappers[idtype] = mapper
                    break
                except KeyError:
                    pass
            if mapper is None:
                raise UnknownIdType("{} is unknown".format(str(idtype)))
        return mapper


class BestIdPerEntryOnlyMixin:
    def filter_best_id(self, arr):
        df = pd.DataFrame(arr)
        order = self.canonical_source_order()
        df["ord"] = df["XRefSource"].apply(lambda x: order.index(self.xrefEnum(x)))
        df = df.sort_values(by=["EntryNr", "ord"]).drop_duplicates(
            subset="EntryNr", keep="first"
        )
        return df.drop(columns="ord").to_records(
            index=False, column_dtypes=dict(arr.dtype.descr)
        )


class DescriptionSearcher(object):
    def __init__(self, db):
        self.entry_tab = db.get_hdf5_handle().get_node("/Protein/Entries")
        self.desc_index = SuffixSearcher.from_tablecolumn(
            self.entry_tab, "DescriptionOffset"
        )

    @timethis(logging.DEBUG)
    def search_term(self, term, limit=None):
        return self.desc_index.find(term, limit=limit)


class FastMapper(object):
    """GO Function projection to sequences from OMA hdf5 file"""

    def __init__(self, db):
        self.db = db

    def iter_projected_goannotations(self, records, target_species=None):
        # gene ontology fast mapping, uses exact / approximate search.
        # todo: implement taxonomic restriction.
        # Input: iterable of biopython SeqRecords

        entrynr_range = None
        if target_species is not None:
            entrynr_range = self.db.id_mapper["OMA"].genome_range(target_species)

        for rec in records:
            logger.debug("projecting function to {}".format(rec))
            if len(rec) < 25:
                logger.warning("Skipping short sequence (len={} AA)".format(len(rec)))
                continue
            t0 = time.time()
            r = self.db.seq_search.search(str(rec.seq), entrynr_range=entrynr_range)
            logger.info(
                "sequence matching of {} ({} AA) took {:.3f}sec".format(
                    rec.id, len(rec), time.time() - t0
                )
            )
            if r is not None:
                logger.debug(str(r))
                go_df = None
                if r[0] == "exact":
                    tdfs1 = []
                    for enum in r[1]:
                        df = self.db.get_gene_ontology_annotations(
                            enum, as_dataframe=True
                        )
                        if df is not None:
                            df["With"] = "Exact:{}".format(
                                self.db.id_mapper["Oma"].map_entry_nr(enum)
                            )
                            tdfs1.append(df)
                    if len(tdfs1) > 0:
                        go_df = pd.concat(tdfs1, ignore_index=True)

                else:
                    # Take best match, that has go-annotations unless it's alignment score
                    # is less than 80% of the highest scoring match.
                    max_score = r[1][0][1]["score"]
                    for enr, match_res in r[1]:
                        score = match_res["score"]
                        if score < 70 or score < 0.8 * max_score:
                            break

                        go_df = self.db.get_gene_ontology_annotations(
                            enr, as_dataframe=True
                        )
                        cnt = len(go_df) if go_df is not None else 0
                        logger.debug(
                            "match: enr: {}, score:{}, go_anno: {}".format(
                                enr, score, cnt
                            )
                        )

                        if go_df is not None:
                            go_df["With"] = "Approx:{}:{}".format(
                                self.db.id_mapper["Oma"].map_entry_nr(enr),
                                score,
                            )
                            break
                if go_df is not None:
                    go_df["DB"] = "OMA_FastMap"
                    go_df["Assigned_By"] = go_df["DB"]
                    go_df["DB_Object_ID"] = rec.id
                    go_df["DB_Object_Symbol"] = go_df["DB_Object_ID"]
                    go_df["Evidence"] = "IEA"
                    go_df["DB:Reference"] = "OMA_Fun:002"
                    go_df["Taxon_ID"] = "taxon:-1"
                    len_with_dupl = len(go_df)
                    go_df.drop_duplicates(inplace=True)
                    logger.debug(
                        "cleaning duplicates: from {} to {} annotations".format(
                            len_with_dupl, len(go_df)
                        )
                    )
                    for row in go_df.to_dict("records"):
                        yield row

    def write_annotations(self, file, seqrecords, target_species=None):
        """Project annotations and write them to file

        This method takes a filehandle and an iterable of BioPython
        SeqRecords objects as input. The function computes the
        projected annotations and writes them to the file in gaf
        format.

        :param file: filehandle to write annotations to
        :param seqrecords: input sequencs to project functions to
        """

        file.write("!gaf-version: {}\n".format(GAF_VERSION))
        file.write("!Project Name: OMA Fast Function Projection\n")
        file.write("!Date created: {}\n".format(time.strftime("%c")))
        file.write("!Contact Email: contact@omabrowser.org\n")
        for anno in self.iter_projected_goannotations(seqrecords, target_species):
            GOA.writerec(anno, file, GOA.GAF20FIELDS)


HogMapperTuple = collections.namedtuple(
    "HogMapperTuple", "query, closest_entry_nr, target, distance, score"
)


class ClosestSeqMapper(object):
    """Fast mapper of sequeces to closest sequence and taken their HOG annotation"""

    def __init__(self, db, threaded=False):
        self.db = Database(db.get_hdf5_handle().filename) if threaded else db

    def _get_main_protein_from_entrynr(self, enr):
        entry = ProteinEntry(self.db, self.db.entry_by_entry_nr(enr))
        return entry.get_main_isoform()

    def imap_sequences(self, seqrecords, target_species=None):
        """maps an iterator of Bio.Seq.SeqRecords to database and
        yields the mappings."""
        entrynr_range = None
        if target_species is not None:
            entrynr_range = self.db.id_mapper["OMA"].genome_range(target_species)

        for rec in seqrecords:
            logger.debug("projecting {} to closest sequence".format(rec))
            r = self.db.seq_search.search(
                str(rec.seq), compute_distance=True, entrynr_range=entrynr_range
            )
            if r is not None:
                if r[0] == "exact":
                    # r[1] contains list of entry nr with exact sequence match (full length)
                    for enr in r[1]:
                        p = self._get_main_protein_from_entrynr(enr)
                        # TODO: good score for identical sequences
                        yield HogMapperTuple(rec.id, enr, p, 0, 0)
                else:
                    candidate_matches = r[1]
                    if len(candidate_matches) <= 0:
                        continue
                    # Take best matches, up to score>80 & score_i > .5*score_max
                    score_max = candidate_matches[0][1]["score"]
                    for enr, match_res in candidate_matches:
                        if (
                            match_res["score"] < 80
                            or match_res["score"] < 0.5 * score_max
                        ):
                            break
                        p = self._get_main_protein_from_entrynr(enr)
                        # score, distance, distvar = pyopa.MutipleAlEnv
                        yield HogMapperTuple(
                            rec.id, enr, p, match_res["distance"], match_res["score"]
                        )
