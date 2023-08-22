from __future__ import division, print_function

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
import resource
import subprocess
import time
from builtins import str, chr, range, object, super, bytes
from tempfile import NamedTemporaryFile

# from pebble import ProcessExpired, ProcessPool

import familyanalyzer
import lxml.html
import numpy
import numpy.lib.recfunctions
import pandas
import tables
from PySAIS import sais
from future.standard_library import hooks
from tqdm import tqdm

from . import locus_parser
from . import suffixsearch
from . import tablefmt
from .KmerEncoder import KmerEncoder
from .OrthoXMLSplitter import OrthoXMLSplitter
from .geneontology import GeneOntology, OntologyParser, FreqAwareGeneOntology
from .homoeologs import HomeologsConfidenceCalculator
from .synteny import SyntenyScorer
from . import hoghelper
from .. import common, version

with hooks():
    import urllib.request

hog_re = re.compile(
    rb"((?P<prefix>HOG):)?(?P<rel>[A-Z]?)(?P<fam>\d+)(\.(?P<subfamid>[0-9a-z.]+))?$"
)


class DarwinException(Exception):
    def __init__(self, stderr, stdout):
        msg = ["Exception running darwin. End of stdout:\n"]
        if len(stdout) > 100:
            msg.append(" ...[{:d}]...".format(len(stdout) - 100))
        msg.append("{!r}".format(stdout[-100:]))
        msg.append("End of stderr:\n")
        if len(stderr) > 100:
            msg.append(" ...[{:d}]...".format(len(stderr) - 200))
        msg.append("{!r}".format(stderr[-200:]))
        super(DarwinException, self).__init__("".join(msg))
        self.stderr = stderr
        self.stdout = stdout


def callDarwinExport(func, drwfile=None):
    """Function starts a darwin session, loads convert.drw file
    and calls the darwin function passed as argument. The output
    is expected to be written by darwin in json format into the
    file specified by 'outfn'.
    This function returns the parsed json datastructure"""

    with NamedTemporaryFile(suffix=".dat") as tmpfile:
        if drwfile is None:
            drwfile = os.path.abspath(os.path.splitext(__file__)[0] + ".drw")
        # with open(os.devnull, 'w') as DEVNULL:
        stacksize = resource.getrlimit(resource.RLIMIT_STACK)
        common.package_logger.debug("current stacklimit: {}".format(stacksize))
        common.package_logger.debug(
            "setting stacklimit: {}".format((max(stacksize) - 1, stacksize[1]))
        )
        resource.setrlimit(resource.RLIMIT_STACK, (min(stacksize), stacksize[1]))
        p = subprocess.Popen(
            ["darwin", "-q", "-E"],
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        drw_cmd = "outfn := '{}': ReadProgram('{}'): {}; done;".format(
            tmpfile.name, drwfile, func
        ).encode("utf-8")
        common.package_logger.debug("calling darwin function: {}".format(func))
        stdout, stderr = p.communicate(input=drw_cmd)
        if p.returncode > 0:
            raise DarwinException(stderr, stdout)

        trans_tab = "".join(str(chr(x)) for x in range(128)) + " " * 128
        if not os.path.exists(tmpfile.name) or os.path.getsize(tmpfile.name) == 0:
            raise DarwinException(stderr, stdout)
        with open(tmpfile.name, "r") as jsonData:
            rawdata = jsonData.read()
            return json.loads(rawdata.translate(trans_tab))


def uniq(seq, transform=None):
    """return uniq elements of a list, preserving order

    The transform argument can be used to check the existance
    on a transformed form of the elements, e.g. only to look
    at a subset of the attributes/values.

    :Example:

        >>> seq = [(1,2,3),(1,2,4),(2,2,2)]
        >>> uniq(seq, transform=lambda x: x[:1])
        [(1, 2, 3), (2, 2, 2)]
        >>> uniq(seq)
        [(1, 2, 3), (1, 2, 4), (2, 2, 2)]


    :param seq: an iterable to be analyzed
    """

    def pass_trough(x):
        return x

    seen = set()
    if transform is None:
        transform = pass_trough
    return [x for x in seq if not (transform(x) in seen or seen.add(transform(x)))]


def silentremove(filename):
    """Function to remove a given file. No exception is raised if the
    file does not exist. Other errors are passed to the user.
    :param filename: the path of the file to be removed"""
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured


def gz_is_empty(fname):
    """Test if gzip file fname is empty

    Return True if the uncompressed data in fname has zero length
    or if fname itself has zero length
    Raises OSError if fname has non-zero length and is not a gzip file
    """
    with gzip.open(fname, "rb") as f:
        data = f.read(1)
    return len(data) == 0


def load_tsv_to_numpy(args):
    fn, off1, off2, swap = args
    rel_desc = tablefmt.PairwiseRelationTable
    # we need to get the enum as a dict to be able to extend it
    # with the reversed labels, i.e. n:1
    relEnum = rel_desc.columns["RelType"].enum._names
    relEnum["n:1"] = relEnum["m:1"]
    relEnum["1:m"] = relEnum["1:n"]
    relEnum["n:m"] = relEnum["m:n"]
    read_dir = -1 if swap else 1
    tsv_dtype = [
        ("EntryNr1", "i4"),
        ("EntryNr2", "i4"),
        ("Score", "f4"),
        ("RelType", "i1"),
        ("AlignmentOverlap", "f2"),
        ("Distance", "f4"),
    ]
    for curNr, curFn in enumerate([fn, fn.replace(".ext.", ".")]):
        try:
            if gz_is_empty(curFn):
                return numpy.empty(0, dtype=tables.dtype_from_descr(rel_desc))
            with gzip.GzipFile(curFn) as fh:
                data = numpy.genfromtxt(
                    fh,
                    dtype=tsv_dtype,
                    names=[_[0] for _ in tsv_dtype],
                    delimiter="\t",
                    usecols=(0, 1, 2, 3, 4, 5),
                    converters={
                        "EntryNr1": lambda nr: int(nr) + off1,
                        "EntryNr2": lambda nr: int(nr) + off2,
                        "RelType": lambda rel: (
                            relEnum[rel[::read_dir].decode()]
                            if len(rel) <= 3
                            else relEnum[rel.decode()]
                        ),
                        "Score": lambda score: float(score) / 100,
                    },
                )
                break
        except OSError as e:
            if curNr < 1:
                common.package_logger.info("tried to load {}".format(curFn))
                pass
            else:
                raise e

    if swap:
        reversed_cols = tuple(data.dtype.names[z] for z in (1, 0, 2, 3, 4, 5))
        data.dtype.names = reversed_cols
    full_table = numpy.empty(data.size, dtype=tables.dtype_from_descr(rel_desc))
    common_cols = list(data.dtype.names)
    full_table[common_cols] = data[common_cols]
    for col_not_in_tsv in set(full_table.dtype.names) - set(data.dtype.names):
        full_table[col_not_in_tsv] = rel_desc.columns[col_not_in_tsv].dflt
    return full_table


def read_vps_from_tsv(gs, ref_genome, basedir=None):
    ref_genome_idx = gs.get_where_list(
        "(UniProtSpeciesCode==code)", condvars={"code": ref_genome}
    )[0]
    job_args = []
    if basedir is None:
        basedir = os.path.join(os.environ["DARWIN_OMADATA_PATH"], "Phase4")
    for g in range(len(gs)):
        if g == ref_genome_idx:
            continue
        g1, g2 = sorted((g, ref_genome_idx))
        off1, off2 = gs.read_coordinates(numpy.array((g1, g2)), "EntryOff")
        fn = os.path.join(
            basedir,
            gs.cols.UniProtSpeciesCode[g1].decode(),
            gs.cols.UniProtSpeciesCode[g2].decode() + ".orth.txt.gz",
        )
        tup = (fn, off1, off2, g1 != ref_genome_idx)
        common.package_logger.debug("adding job: {}".format(tup))
        job_args.append(tup)

    pool = mp.Pool(processes=min(os.cpu_count(), 12))
    all_pairs = pool.map(load_tsv_to_numpy, job_args)
    pool.close()
    common.package_logger.info("loaded vps for {}".format(ref_genome.decode()))
    return numpy.lib.recfunctions.stack_arrays(all_pairs, usemask=False)


def load_hogs_at_level(fname, level):
    with tables.open_file(fname, "r") as h5:
        lev = level.encode("utf-8") if isinstance(level, str) else level
        tab = h5.get_node("/HogLevel")
        extended_dtype = numpy.dtype(tab.dtype.descr + [("HogLevelRowIdx", "i4")])
        hog_it = (
            tuple(row.fetch_all_fields()) + (row.nrow,)
            for row in tab.where("Level == lev")
        )
        hogs = numpy.fromiter(hog_it, dtype=extended_dtype)
        hogs.sort(order="ID")
        hogs["IdxPerLevelTable"] = numpy.arange(len(hogs))
        return hogs


class DataImportError(Exception):
    pass


def _load_taxonomy_without_ref_to_itselfs(data):
    dtype = tables.dtype_from_descr(tablefmt.TaxonomyTable)
    arr = numpy.array([tuple(x) for x in data], dtype=dtype)
    clean = arr[numpy.where(arr["NCBITaxonId"] != arr["ParentTaxonId"])]
    return clean


def compute_ortholog_types(data, genome_offs):
    """this function computes the type of orthologs from the data and sets in
    the RelType column.

    :param data: a numpy recarray corresponding to the `numpy.dtype` of
           `tablefmt.PairwiseRelationTable`
    :param genome_offs: a numpy array with the genome offsets, i.e. the entry
           numbers where the next genome starts

    :returns: a modified version of data
    """
    typEnum = tablefmt.PairwiseRelationTable.columns.get("RelType").enum
    query_type = {
        val: "m" if cnt > 1 else "1"
        for val, cnt in zip(*numpy.unique(data["EntryNr2"], return_counts=True))
    }

    def genome_idx(enr):
        return numpy.searchsorted(genome_offs, enr - 1, side="right")

    g0 = genome_idx(data[0]["EntryNr2"])
    it = numpy.nditer(data, flags=["c_index"], op_flags=["readwrite"])
    while not it.finished:
        row0 = it[0]
        i1 = it.index + 1
        # we move i1 forward to the row where the next genome starts, i.e. the
        # current query changes the species or the query itself changes
        while i1 < len(data):
            row1 = data[i1]
            g1 = genome_idx(row1["EntryNr2"])
            if g1 != g0 or row0["EntryNr1"] != row1["EntryNr1"]:
                break
            i1 += 1
        subj_type = "n" if i1 - it.index > 1 else "1"
        while not it.finished and it.index < i1:
            typ = "{}:{}".format(query_type[int(it[0]["EntryNr2"])], subj_type)
            it[0]["RelType"] = typEnum[typ]
            it.iternext()
        g0 = g1


def get_or_create_tables_node(h5, path, desc=None):
    """return the node of a given path from the h5 file

    If the node does not yet exist, it is created (including potential
    inexistant internal nodes).

    :param h5: Handle to the hdf5 object
    :param str path: Path of the node to return
    :param str desc: Description to be added to the node"""
    try:
        grp = h5.get_node(path)
    except tables.NoSuchNodeError:
        base, name = os.path.split(path)
        grp = h5.create_group(base, name, title=desc, createparents=True)
    return grp


def create_index_for_columns(tab, *cols):
    if not isinstance(tab, tables.Table):
        raise TypeError("tab argument must be table node")
    for col in cols:
        if not tab.colindexed[col]:
            tab.colinstances[col].create_csindex()
        else:
            tab.colinstances[col].reindex_dirty()


def create_fast_famhoglevel_lookup(hoglevtab):
    max_fam_nr = hoglevtab[hoglevtab.colindexes["Fam"][-1]]["Fam"]
    lookup = numpy.zeros(max_fam_nr + 1, dtype=[("start", "i4"), ("stop", "i4")])
    lst_fam = 0
    lst_start = 0
    for row in hoglevtab.iterrows():
        if row["Fam"] != lst_fam:
            lookup[lst_fam] = (lst_start, row.nrow)
            lst_fam = row["Fam"]
            lst_start = row.nrow
    lookup[lst_fam] = (lst_start, row.nrow + 1)
    return lookup


def create_and_store_fast_famhoglevel_lookup(h5, hoglevtab, array_path):
    try:
        n = h5.get_node(array_path)
        n.remove()
    except tables.NoSuchNodeError:
        pass
    lookup = create_fast_famhoglevel_lookup(hoglevtab)
    node, name = os.path.split(array_path)
    h5.create_table(node, name, expectedrows=len(lookup), obj=lookup)


class DarwinExporter(object):
    DB_SCHEMA_VERSION = "3.6"
    DRW_CONVERT_FILE = os.path.abspath(os.path.splitext(__file__)[0] + ".drw")

    def __init__(self, path, logger=None, mode=None):
        self.logger = logger if logger is not None else common.package_logger
        fn = os.path.normpath(
            os.path.join(os.getenv("DARWIN_BROWSERDATA_PATH", ""), path)
        )
        if mode is None:
            mode = "append" if os.path.exists(fn) else "write"
        self._compr = tables.Filters(complevel=6, complib="zlib", fletcher32=True)
        self.h5 = tables.open_file(fn, mode=mode[0], filters=self._compr)
        self.logger.info(
            "opened {} in {} mode, options {} ; pyoma {}".format(
                fn, mode, str(self._compr), version()
            )
        )
        if mode == "write":
            self.h5.root._f_setattr("convertion_start", time.strftime("%c"))
            self.h5.root._f_setattr("pyoma_version", version())

    def call_darwin_export(self, func):
        return callDarwinExport(func, self.DRW_CONVERT_FILE)

    def _get_or_create_node(self, path, desc=None):
        return get_or_create_tables_node(self.h5, path, desc)

    def create_table_if_needed(self, parent, name, drop_data=False, **kwargs):
        """create a table if needed.

        The function only checks whether a table exists with that name,
        but not if it is compatible with the passed arguments.
        if you pass data with the `obj` argument, this data is appended
        to the table. If you set `drop_data` to True, data that was
        previously in the existing table is dropped prior to adding new
        data."""
        try:
            tab = self.h5.get_node(parent, name=name)
            if drop_data:
                tab.remove_rows(0, tab.nrows)
            if "obj" in kwargs:
                tab.append(kwargs["obj"])
        except tables.NoSuchNodeError:
            tab = self.h5.create_table(parent, name, **kwargs)
        return tab

    def get_version(self):
        """return version of the dataset.

        Default implementation searches for 'mname' in Matrix or matrix_stats.drw files.
        """
        for fname in ("Matrix", "matrix_stats.drw"):
            with open(
                os.path.join(os.environ["DARWIN_BROWSERDATA_PATH"], fname), "r"
            ) as fh:
                for i, line in enumerate(fh):
                    if line.startswith("mname :="):
                        match = re.match(r"mname := \'(?P<version>[^\']*)\'", line)
                        return match.group("version")
                    if i > 1000:
                        break
        raise DataImportError("No version information found")

    def add_version(self, release_char=None):
        if release_char is not None:
            if not re.match(r"^[A-Z]*$", release_char):
                raise ValueError(
                    "Unexpected release_char argument: {}. Must be a capital ascii character".format(
                        release_char
                    )
                )
        else:
            release_char = ""
        version = self.get_version()
        self.h5.set_node_attr("/", "oma_version", version)
        self.h5.set_node_attr("/", "oma_release_char", release_char)
        self.h5.set_node_attr("/", "pytables", tables.get_pytables_version())
        self.h5.set_node_attr("/", "hdf5_version", tables.get_hdf5_version())
        self.h5.set_node_attr("/", "db_schema_version", self.DB_SCHEMA_VERSION)

    def add_species_data(self):
        cache_file = os.path.join(
            os.getenv("DARWIN_NETWORK_SCRATCH_PATH", ""), "pyoma", "gs.json"
        )
        if os.path.exists(cache_file):
            with open(cache_file, "r") as fd:
                data = json.load(fd)
        else:
            data = self.call_darwin_export("GetGenomeData();")
        gstab = self.h5.create_table(
            "/", "Genome", tablefmt.GenomeTable, expectedrows=len(data["GS"])
        )
        gs_data = self._parse_date_columns(data["GS"], gstab)
        self._write_to_table(gstab, gs_data)
        create_index_for_columns(gstab, "NCBITaxonId", "UniProtSpeciesCode", "EntryOff")

        taxtab = self.h5.create_table(
            "/", "Taxonomy", tablefmt.TaxonomyTable, expectedrows=len(data["Tax"])
        )
        self._write_to_table(taxtab, _load_taxonomy_without_ref_to_itselfs(data["Tax"]))
        create_index_for_columns(taxtab, "NCBITaxonId")

    def _parse_date_columns(self, data, tab):
        """convert str values in a date column to epoch timestamps"""
        time_cols = [
            i for i, col in enumerate(tab.colnames) if tab.coldescrs[col].kind == "time"
        ]
        dflts = [tab.coldflts[col] for col in tab.colnames]

        def map_data(col, data):
            try:
                val = data[col]
                if col in time_cols and isinstance(val, str):
                    if val == "":
                        return dflts[col]
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
                return val
            except IndexError:
                return dflts[col]

        arr = numpy.empty(len(data), dtype=tab.dtype)
        for i, row in enumerate(data):
            as_tup = tuple(map_data(c, row) for c in range(len(dflts)))
            arr[i] = as_tup
        return arr

    def _convert_to_numpyarray(self, data, tab):
        """convert a list of list dataset into a numpy rec array that
        corresponds to the table definition of `tab`.

        :param data: the data to be converted.
        :param tab: a pytables table node."""

        enum_cols = {
            i: tab.get_enum(col)
            for (i, col) in enumerate(tab.colnames)
            if tab.coltypes[col] == "enum"
        }
        dflts = [tab.coldflts[col] for col in tab.colnames]

        def map_data(col, data):
            try:
                val = data[col]
                return enum_cols[col][val]
            except IndexError:
                return dflts[col]
            except KeyError:
                return val

        arr = numpy.empty(len(data), dtype=tab.dtype)
        for i, row in enumerate(data):
            as_tup = tuple(map_data(c, row) for c in range(len(dflts)))
            arr[i] = as_tup
        return arr

    def add_orthologs(self):
        genome_offs = self.h5.root.Genome.col("EntryOff")
        anygenome = self.h5.root.Genome[0]["UniProtSpeciesCode"].decode()
        basedir = None
        for base in ("DARWIN_OMADATA_PATH", "DARWIN_OMA_SCRATCH_PATH"):
            testdir = os.path.join(os.getenv(base, "/"), "Phase4", anygenome)
            if os.path.isdir(testdir) and any(
                map(lambda x: x.endswith(".orth.txt.gz"), os.listdir(testdir))
            ):
                basedir = os.path.join(os.getenv(base), "Phase4")
                break

        self.logger.info("using {} as base dir for pairwise orthology".format(basedir))
        for gs in self.h5.root.Genome.iterrows():
            genome = gs["UniProtSpeciesCode"].decode()
            rel_node_for_genome = self._get_or_create_node(
                "/PairwiseRelation/{}".format(genome)
            )
            if "VPairs" not in rel_node_for_genome:
                cache_file = os.path.join(
                    os.getenv("DARWIN_NETWORK_SCRATCH_PATH", ""),
                    "pyoma",
                    "vps",
                    "{}.json".format(genome),
                )
                if os.path.exists(cache_file):
                    with open(cache_file, "r") as fd:
                        data = json.load(fd)
                elif basedir is not None:
                    # we have the *.orth.txt.gz files in basedir
                    data = read_vps_from_tsv(
                        self.h5.root.Genome, genome.encode("utf-8"), basedir=basedir
                    )
                else:
                    # fallback to read from VPsDB
                    data = self.call_darwin_export("GetVPsForGenome({})".format(genome))

                vp_tab = self.h5.create_table(
                    rel_node_for_genome,
                    "VPairs",
                    tablefmt.PairwiseRelationTable,
                    expectedrows=len(data),
                )
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, vp_tab)
                if numpy.any(
                    data["RelType"]
                    >= tablefmt.PairwiseRelationTable.columns.get("RelType").enum["n/a"]
                ):
                    compute_ortholog_types(data, genome_offs)
                self._write_to_table(vp_tab, data)
                create_index_for_columns(vp_tab, "EntryNr1")

    def add_same_species_relations(self):
        for gs in self.h5.root.Genome.iterrows():
            genome = gs["UniProtSpeciesCode"].decode()
            rel_node_for_genome = self._get_or_create_node(
                "/PairwiseRelation/{}".format(genome)
            )
            if "within" not in rel_node_for_genome:
                cache_file = os.path.join(
                    os.getenv("DARWIN_NETWORK_SCRATCH_PATH", ""),
                    "pyoma",
                    "cps",
                    "{}.json".format(genome),
                )
                if os.path.exists(cache_file):
                    with open(cache_file, "r") as fd:
                        data = json.load(fd)
                else:
                    # fallback to read from VPsDB
                    data = self.call_darwin_export(
                        "GetSameSpeciesRelations({})".format(genome)
                    )

                ss_tab = self.h5.create_table(
                    rel_node_for_genome,
                    "within",
                    tablefmt.PairwiseRelationTable,
                    expectedrows=len(data),
                )
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, ss_tab)
                self._write_to_table(ss_tab, data)
                create_index_for_columns(ss_tab, "EntryNr1")

    def add_synteny_scores(self):
        """add synteny scores of pairwise relations to database.

        Current implementation only computes synteny scores for
        homoeologs, but easy to extend. Question is rather if we
        need synteny scores for all genome pairs, and if not, how
        to select.

        The computations of the scores are done using :mod:`synteny`
        module of this package."""
        # TODO: compute for non-homoeologs relation as well.
        self.logger.info("Adding synteny scores for polyploid genomes")
        polyploid_genomes = self.h5.root.Genome.where("IsPolyploid==True")
        for genome in polyploid_genomes:
            genome_code = genome["UniProtSpeciesCode"].decode()
            self.logger.info("compute synteny score for {}".format(genome_code))
            synteny_scorer = SyntenyScorer(self.h5, genome_code)
            rels = synteny_scorer.compute_scores()
            self._callback_store_rel_data(
                genome_code, rels, [("SyntenyConservationLocal", "mean_synteny_score")]
            )

    def add_homoeology_confidence(self):
        """adds the homoeology confidence scores to the database.

        This method should be called only after the synteny scores have
        been computed and added to the database.

        The computations are done using :mod:`homoeologs` module."""
        self.logger.info("Adding homoeolog confidence scores")
        polyploid_genomes = self.h5.root.Genome.where("IsPolyploid==True")
        for genome in polyploid_genomes:
            genome_code = genome["UniProtSpeciesCode"].decode()
            self.logger.info("compute homoeolog confidence for {}".format(genome_code))
            homoeolg_scorer = HomeologsConfidenceCalculator(self.h5, genome_code)
            rels = homoeolg_scorer.calculate_scores()
            self._callback_store_rel_data(
                genome_code, rels, [("Confidence", "fuzzy_confidence_scaled")]
            )

    def _callback_store_rel_data(self, genome, rels_df, assignments):
        tab = self.h5.get_node("/PairwiseRelation/{}/within".format(genome))
        df_all = pandas.DataFrame(tab.read())
        if "entry_nr1" in list(rels_df):
            enr_col_names = ["entry_nr1", "entry_nr2"]
        else:
            enr_col_names = ["EntryNr1", "EntryNr2"]
        merged = pandas.merge(
            df_all,
            rels_df,
            how="left",
            left_on=["EntryNr1", "EntryNr2"],
            right_on=enr_col_names,
            validate="one_to_one",
        )

        for target, source in assignments:
            # replace NaN in column from rels_df by the default value of the target column
            merged.loc[merged[source].isnull(), source] = tab.coldescrs[target].dflt
            # update the data in the target hdf5 column by the source column data
            tab.modify_column(column=merged[source].to_numpy(), colname=target)
        tab.flush()

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
        protGrp = self._get_or_create_node(
            "/Protein", "Root node for protein (oma entries) information"
        )
        protTab = self.h5.create_table(
            protGrp, "Entries", tablefmt.ProteinTable, expectedrows=nrProt
        )
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
                data = self.call_darwin_export(
                    "GetProteinsForGenome({})".format(genome)
                )

            if len(data["seqs"]) != gs["TotEntries"]:
                raise DataImportError(
                    "number of entries ({:d}) does "
                    "not match number of seqs ({:d}) for {}".format(
                        len(data["seqs"]), gs["TotEntries"], genome
                    )
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

                seqOff += self._add_sequence(
                    data["seqs"][nr], protTab.row, seqArr, seqOff
                )
                cdnaOff += self._add_sequence(
                    data["cdna"][nr], protTab.row, cdnaArr, cdnaOff, "CDNA"
                )

                protTab.row["Chromosome"] = data["chrs"][nr]
                protTab.row["AltSpliceVariant"] = data["alts"][nr]
                protTab.row["OmaHOG"] = b" "  # will be assigned later
                protTab.row["CanonicalId"] = b" "  # will be assigned later
                if (
                    protTab.row["AltSpliceVariant"] == 0
                    or protTab.row["AltSpliceVariant"] == protTab.row["EntryNr"]
                ):
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
                self.logger.warning(
                    "{} missmatches in exon-lengths compared to locus info".format(
                        cnt_missmatch_locus
                    )
                )
            for n in (protTab, seqArr, locTab):
                if n.size_in_memory != 0:
                    self.logger.info(
                        "worte %s: compression ratio %3f%%"
                        % (n._v_pathname, 100 * n.size_on_disk / n.size_in_memory)
                    )
        create_index_for_columns(protTab, "EntryNr", "MD5ProteinHash")

    def _write_to_table(self, tab, data):
        if len(data) > 0:
            tab.append(data)
            self.logger.info(
                "wrote %s : compression ratio %.3f%%"
                % (tab._v_pathname, 100 * tab.size_on_disk / tab.size_in_memory)
            )
        else:
            self.logger.info("no data written for {}".format(tab._v_pathname))

    def add_hogs(self, hog_path=None, hog_file=None, tree_filename=None):
        """adds the HOGs to the database

        :param str hog_path: optional, directory where the split HOG files are stored or
                             should be stored. If directory does not exist, the hog_file
                             is split automatically into individual files and stored there.

        :param str hog_file: File containing all HOGs. if hog_path does not
                             exist, this file is split and stored in hog_path.

        :param str tree_filename: newick species tree file."""
        if hog_path is None:
            hog_path = os.path.normpath(
                os.path.join(
                    os.environ["DARWIN_NETWORK_SCRATCH_PATH"], "pyoma", "split_hogs"
                )
            )
        entryTab = self.h5.get_node("/Protein/Entries")
        try:
            release = self.h5.get_node_attr("/", "oma_release_char")
        except AttributeError:
            release = ""
        if tree_filename is None:
            tree_filename = os.path.join(
                os.environ["DARWIN_BROWSERDATA_PATH"], "speciestree.nwk"
            )
        if not os.path.exists(hog_path):
            if hog_file is None:
                hog_file = os.path.join(
                    os.environ["DARWIN_BROWSERDATA_PATH"],
                    "..",
                    "downloads",
                    "oma-hogs.orthoXML.gz",
                )
            splitter = OrthoXMLSplitter(
                hog_file, cache_dir=hog_path, release_char=release
            )
            splitter()
        tax_tab = self.h5.get_node("/Taxonomy")
        tax_2_code = {
            int(row["NCBITaxonId"]): row["UniProtSpeciesCode"].decode()
            for row in self.h5.get_node("/Genome")
        }
        hog_converter = HogConverter(entryTab, release, tax_tab, tax_2_code)
        hog_converter.attach_newick_taxonomy(tree_filename)
        hogTab = self.h5.create_table(
            "/",
            "HogLevel",
            tablefmt.HOGsTable,
            "nesting structure for each HOG",
            expectedrows=1e8,
        )
        self.orthoxml_buffer = self.h5.create_earray(
            "/OrthoXML",
            "Buffer",
            tables.StringAtom(1),
            (0,),
            "concatenated orthoxml files",
            expectedrows=1e9,
            createparents=True,
        )
        self.orthoxml_buffer_augmented = self.h5.create_earray(
            "/OrthoXML",
            "BufferAugmented",
            tables.StringAtom(1),
            (0,),
            "concatenated augmented orthoxml files",
            expectedrows=1e9,
            createparents=True,
        )
        self.orthoxml_index = self.h5.create_table(
            "/OrthoXML",
            "Index",
            tablefmt.OrthoXmlHogTable,
            "Range index per HOG into OrthoXML Buffer",
            expectedrows=5e6,
        )
        for root, dirs, filenames in os.walk(hog_path):
            for fn in filenames:
                if fn.endswith(".augmented"):
                    continue
                try:
                    input_file = os.path.join(root, fn)
                    out_orthoxml = input_file + ".augmented"
                    levels = hog_converter.convert_file(input_file, store=out_orthoxml)
                    hogTab.append(levels)
                    fam_nrs = set([z[0] for z in levels])
                    self.add_orthoxml(input_file, out_orthoxml, fam_nrs)
                except Exception as e:
                    self.logger.error("an error occured while processing " + fn + ":")
                    self.logger.exception(e)
        # flushing index table
        self.orthoxml_index.flush()
        hog_converter.write_hogs()

    def add_orthoxml(self, orthoxml_path, augmented_orthoxml_path, fam_nrs):
        """append orthoxml file content to orthoxml_buffer array and add index for the HOG family"""
        if len(fam_nrs) > 1:
            self.logger.warning(
                "expected only one family per HOG file, but found {}: {}".format(
                    len(fam_nrs), fam_nrs
                )
            )
            self.logger.warning(
                " --> the orthoxml files per family will be not correct, "
                "i.e. they will contain all families of this file."
            )
        offset = {}
        length = {}
        for typ, fpath in zip(
            ("base", "augmented"), (orthoxml_path, augmented_orthoxml_path)
        ):
            buffer = (
                self.orthoxml_buffer
                if typ == "base"
                else self.orthoxml_buffer_augmented
            )
            offset[typ] = len(buffer)
            try:
                with open(fpath, "r") as fh:
                    orthoxml = fh.read().encode("utf-8")
                    length[typ] = len(orthoxml)
                    buffer.append(
                        numpy.ndarray(
                            (length[typ],), buffer=orthoxml, dtype=tables.StringAtom(1)
                        )
                    )
            except IOError:
                length[typ] = 0
        for fam in fam_nrs:
            row = self.orthoxml_index.row
            row["Fam"] = fam
            row["HogBufferOffset"] = offset["base"]
            row["HogBufferLength"] = length["base"]
            row["HogAugmentedBufferOffset"] = offset["augmented"]
            row["HogAugmentedBufferLength"] = length["augmented"]
            row.append()

    def add_cache_of_hogs_by_level(self, nr_procs=None):
        self.logger.info("createing cached HogLevel table per level")
        hl_tab = self.h5.get_node("/HogLevel")
        temp_hoglevel_file = os.path.join(
            os.getenv("DARWIN_NETWORK_SCRATCH_PATH"), "tmp-hoglevel.h5"
        )
        with tables.open_file(temp_hoglevel_file, "w") as hlfh:
            hl_tab._f_copy(hlfh.root)
            create_index_for_columns(hlfh.get_node("/HogLevel"), "Level")

        rel_levels = set(hl_tab.read(field="Level"))
        self.logger.info(
            "found {} levels, start extracting hogs in parallel".format(len(rel_levels))
        )
        lev2tax = {
            row["Name"]: int(row["NCBITaxonId"])
            for row in self.h5.get_node("/Taxonomy").read()
        }
        lev2tax[b"LUCA"] = 0
        idx_per_level = numpy.zeros(len(hl_tab), "i4")
        with concurrent.futures.ProcessPoolExecutor(max_workers=nr_procs) as pool:
            future_to_level = {
                pool.submit(
                    load_hogs_at_level, fname=temp_hoglevel_file, level=level
                ): level
                for level in rel_levels
            }
            for future in concurrent.futures.as_completed(future_to_level):
                level = future_to_level[future]
                try:
                    hogs = future.result()
                    # fallback to level if taxid is not known
                    tab_name = "tax{}".format(lev2tax.get(level, level.decode()))
                    tab = self.h5.create_table(
                        where=f"/AncestralGenomes/{tab_name}",
                        name="Hogs",
                        title="cached HogLevel data for {}".format(level.decode()),
                        obj=hogs,
                        createparents=True,
                        expectedrows=len(hogs),
                    )
                    create_index_for_columns(
                        tab, "Fam", "IsRoot", "NrMemberGenes", "CompletenessScore"
                    )
                    idx_per_level[hogs["HogLevelRowIdx"]] = hogs["IdxPerLevelTable"]
                except Exception as exc:
                    msg = "cannot store cached hogs for {}".format(level)
                    self.logger.exception(msg)
                    pass
        hl_tab.modify_column(column=idx_per_level, colname="IdxPerLevelTable")

    def xref_databases(self):
        return os.path.join(os.environ["DARWIN_BROWSERDATA_PATH"], "ServerIndexed.db")

    def add_xrefs(self):
        self.logger.info("start extracting XRefs, EC and GO annotations")
        db_parser = DarwinDbEntryParser()
        xref_tab = self.h5.create_table(
            "/",
            "XRef",
            tablefmt.XRefTable,
            "Cross-references of proteins to external ids / descriptions",
            expectedrows=1e8,
        )

        ec_tab = self.h5.create_table(
            "/Annotations",
            "EC",
            tablefmt.ECTable,
            "Enzyme Commission annotations",
            expectedrows=1e7,
            createparents=True,
        )
        gs = self.h5.get_node("/Genome").read()
        approx_adder = ApproximateXRefImporter(
            os.path.join(
                os.getenv("DARWIN_NETWORK_SCRATCH_PATH"), "approximate_xrefs.tsv"
            )
        )
        up_mapped_xrefs = UniProtAdditionalXRefImporter(
            *(
                os.path.join(os.getenv("DARWIN_NETWORK_SCRATCH_PATH"), z)
                for z in ("uniprot_xrefs.tsv", "swissprot_xrefs.tsv")
            )
        )
        with DescriptionManager(
            self.h5, "/Protein/Entries", "/Protein/DescriptionBuffer"
        ) as de_man, GeneOntologyManager(
            self.h5, "/Annotations/GeneOntology", "/Ontologies/GO"
        ) as go_man:
            xref_importer = XRefImporter(
                db_parser,
                gs,
                xref_tab,
                ec_tab,
                go_man,
                de_man,
                approx_adder,
                up_mapped_xrefs,
            )
            files = self.xref_databases()
            dbs_iter = fileinput.input(files=files)
            db_parser.parse_entrytags(dbs_iter)
            xref_importer.flush_buffers()

        # create /XRef_EntryNr_offset array
        enr_col = xref_tab.col("EntryNr")
        if not (enr_col[:-1] <= enr_col[1:]).all():
            raise DataImportError("EntryNr column is not sorted in /XRef")
        enrs = numpy.arange(enr_col[-1] + 2)
        idx = numpy.searchsorted(enr_col, enrs).astype("i4")
        self.h5.create_carray(
            "/",
            "XRef_EntryNr_offset",
            obj=idx,
            title="Array containing the row index in the XRef table at entry_nr: XRef[a[enr]] = first row of EntryNr",
        )

    def add_group_metadata(self):
        m = OmaGroupMetadataLoader(self.h5)
        m.add_data()

    def add_roothog_metadata(self):
        m = RootHOGMetaDataLoader(self.h5)
        m.add_data()

    def close(self):
        self.h5.root._f_setattr("conversion_end", time.strftime("%c"))
        self.h5.flush()
        self.h5.close()
        self.logger.info("closed {}".format(self.h5.filename))

    def create_indexes(self):
        self.logger.info("creating indexes for HogLevel table")
        hogTab = self.h5.get_node("/HogLevel")
        create_index_for_columns(
            hogTab, "Fam", "ID", "Level", "NrMemberGenes", "CompletenessScore", "IsRoot"
        )
        create_and_store_fast_famhoglevel_lookup(
            self.h5, hogTab, "/HogLevel_fam_lookup"
        )

        orthoxmlTab = self.h5.get_node("/OrthoXML/Index")
        create_index_for_columns(orthoxmlTab, "Fam")

        self.logger.info("creating missing indexes for Entries table")
        entryTab = self.h5.get_node("/Protein/Entries")
        create_index_for_columns(
            entryTab, "EntryNr", "OmaHOG", "OmaGroup", "MD5ProteinHash"
        )

        self.logger.info("creating suffix index for Descriptions")
        desc_buffer = self.h5.get_node("/Protein/DescriptionBuffer")
        suffixsearch.create_suffix_index(entryTab, "DescriptionOffset", desc_buffer)

        self.logger.info("creating index for xrefs (EntryNr and XRefId)")
        xrefTab = self.h5.get_node("/XRef")
        create_index_for_columns(xrefTab, "EntryNr", "XRefId")
        self.logger.info("creating suffix index for XRefId")
        suffixsearch.create_suffix_index(xrefTab, "XRefId")

        self.logger.info("creating index for go (EntryNr and TermNr)")
        goTab = self.h5.get_node("/Annotations/GeneOntology")
        create_index_for_columns(goTab, "EntryNr", "TermNr")

        self.logger.info("creating index for EC (EntryNr)")
        ec_tab = self.h5.get_node("/Annotations/EC")
        create_index_for_columns(ec_tab, "EntryNr", "ECacc")

        self.logger.info("creating index for domains (EntryNr)")
        domtab = self.h5.get_node("/Annotations/Domains")
        create_index_for_columns(domtab, "EntryNr", "DomainId")

        self.logger.info(
            "creating indexes for HOG to prevalent domains " "(Fam and DomainId)"
        )
        dom2hog_tab = self.h5.get_node("/HOGAnnotations/Domains")
        create_index_for_columns(dom2hog_tab, "DomainId")
        domprev_tab = self.h5.get_node("/HOGAnnotations/DomainArchPrevalence")
        create_index_for_columns(domprev_tab, "Fam")

        self.logger.info("createing indexes for Domain Descriptions")
        domdesc = self.h5.get_node("/Annotations/DomainDescription")
        create_index_for_columns(domdesc, "DomainId")
        suffixsearch.create_suffix_index(domdesc, "Description")

    def _iter_canonical_xref(self):
        """extract one canonical xref id for each protein.

        We take the first valid xref per gene with the ordering of xrefsources
        as given in the xrefsource_order."""
        xrefsource_order = (
            "UniProtKB/SwissProt",
            "UniProtKB/TrEMBL",
            "Ensembl Gene",
            "Ensembl Protein",
            "FlyBase",
            "WormBase",
            "EnsemblGenomes",
            "RefSeq",
            "SourceID",
        )

        xrefs = self.h5.get_node("/XRef")
        source_enum = xrefs.get_enum("XRefSource")
        canonical_sources = [source_enum[z] for z in xrefsource_order]
        max_acceptable_verif_value = xrefs.get_enum("Verification")["unchecked"]
        current_protein = None
        past_proteins = set([])
        for xref in xrefs:
            if xref["Verification"] > max_acceptable_verif_value:
                continue
            if xref["EntryNr"] != current_protein:
                if current_protein:
                    past_proteins.add(current_protein)
                    yield (current_protein, current_xref[1])
                current_protein = xref["EntryNr"]
                current_xref = (1000, b"")  # init with a sentinel
                if current_protein in past_proteins:
                    raise DataImportError("Data in /XRef is not grouped w.r.t. EntryNr")
            try:
                rank = canonical_sources.index(xref["XRefSource"])
                if rank < current_xref[0]:
                    current_xref = (rank, xref["XRefId"])
            except ValueError:
                pass
        if current_protein:
            yield (current_protein, current_xref[1])

    def add_canonical_id(self):
        """add one canonical xref id to the /Protein/Entries table."""
        self.logger.info("adding canonical ids for each protein...")
        prot_tab = self.h5.get_node("/Protein/Entries")
        canonical_ids = numpy.chararray(
            shape=(len(prot_tab),), itemsize=prot_tab.cols.CanonicalId.dtype.itemsize
        )
        for eNr, canonical_id in self._iter_canonical_xref():
            row_nr = eNr - 1
            row = prot_tab[row_nr]
            if row["EntryNr"] != eNr:
                self.logger.warn(
                    "Entries table not properly sorted: {}, expected {}".format(
                        row["EntryNr"], eNr
                    )
                )
                raise DataImportError("Entries table not properly sorted")
            canonical_ids[row_nr] = canonical_id
        prot_tab.modify_column(
            0, len(prot_tab), 1, column=canonical_ids, colname="CanonicalId"
        )
        prot_tab.flush()

    def add_domain_info(self, domains):
        self.logger.info("adding domain information...")
        domtab = self.h5.create_table(
            "/Annotations",
            "Domains",
            tablefmt.DomainTable,
            createparents=True,
            expectedrows=1e7,
        )
        entrytab = self.h5.get_node("/Protein/Entries")
        md5_to_enr = collections.defaultdict(list)
        for e in entrytab:
            md5_to_enr[e["MD5ProteinHash"]].append(e["EntryNr"])

        buffer = []
        for i, domain in enumerate(domains):
            for entry_nr in md5_to_enr[domain.md5.encode("utf-8")]:
                buffer.append((entry_nr, domain.id, domain.coords))
                if len(buffer) > 5000:
                    domtab.append(buffer)
                    buffer = []
            if i % 50000 == 0:
                self.logger.info("processed {:d} domain annotations so far".format(i))
        if len(buffer) > 0:
            domtab.append(buffer)
        domtab.flush()

    def add_domainname_info(self, domainname_infos):
        self.logger.info("adding domain name information...")
        dom_name_tab = self.h5.create_table(
            "/Annotations",
            "DomainDescription",
            tablefmt.DomainDescriptionTable,
            createparents=True,
            expectedrows=2e5,
        )
        buffer = []
        for i, dom_info in enumerate(domainname_infos):
            buffer.append(dom_info)
            if len(buffer) > 5000:
                self._write_to_table(dom_name_tab, buffer)
                buffer = []
            if i % 50000 == 0:
                self.logger.info(
                    "processed {:d} domain name descriptions so far".format(i)
                )
        if len(buffer) > 0:
            self._write_to_table(dom_name_tab, buffer)
        dom_name_tab.flush()

    def update_summary_stats(self):
        """update the summary statistics of xrefs & go.

        The function analyses the well-known xref sources as well as
        GO annotations and computes aggregated counts for
        all / in OMA Group / in HOGs for all of them.
        """
        for tab_name, sum_fun in [
            ("/Annotations/GeneOntology", self.count_gene_ontology_summary),
            ("/XRef", self.count_xref_summary),
        ]:
            summary = sum_fun()
            tab = self.h5.get_node(tab_name)
            for attr, val in summary.items():
                tab.set_attr(attr, val)

        group_sizes = self.collect_group_sizes()
        summary = self._get_or_create_node("/Summary", "Various Summary Statistics")
        for group_type in group_sizes.keys():
            grp_size_tab = self.create_table_if_needed(
                summary,
                "{}_size_hist".format(group_type),
                description=tablefmt.GroupsizeHistogram,
                drop_data=True,
            )
            data = sorted(group_sizes[group_type].items())
            grp_size_tab.append(data)

        cov_fracs = self.add_domain_covered_sites_counts()
        cov_hist, bins = numpy.histogram(
            cov_fracs[cov_fracs > 0], bins=numpy.linspace(0, 1, 51)
        )
        cov_hist_data = numpy.zeros(50, dtype=[("BinEndValue", "f4"), ("Counts", "i4")])
        cov_hist_data["BinEndValue"] = bins[1:]
        cov_hist_data["Counts"] = cov_hist
        dom_cov_hist_tab = self.create_table_if_needed(
            summary, "Domain_coverage_hist", drop_data=True, obj=cov_hist_data
        )
        dom_cov_hist_tab.set_attr(
            "frac_genes_w_domain", len(cov_fracs[cov_fracs > 0]) / len(cov_fracs)
        )
        dom_cov_hist_tab.set_attr("mean_coverage_overall", numpy.mean(cov_fracs))
        dom_cov_hist_tab.set_attr(
            "mean_coverage_w_domain", numpy.mean(cov_fracs[cov_fracs > 0])
        )

    def collect_goterm_freqs(self):
        self.logger.info("Computing gene ontology annotations term frequencies")
        go_tab = self.h5.get_node("/Annotations/GeneOntology")
        fp = io.StringIO(self.h5.get_node("/Ontologies/GO").read().tobytes().decode())
        gof = FreqAwareGeneOntology(OntologyParser(fp), rels=None)
        gof.parse()

        def iter_entry_terms(go_tab, filter=None):
            for (enr, term), row_iter in itertools.groupby(
                go_tab, operator.itemgetter("EntryNr", "TermNr")
            ):
                if filter is not None:
                    evid = {row["Evidence"] for row in row_iter}
                    if evid not in filter:
                        continue
                yield {"TermNr": int(term)}

        gof.estimate_freqs(iter_entry_terms(go_tab))
        return gof.cnts, gof.tot_cnts

    def count_gene_ontology_summary(self):
        self.logger.info("Bulding gene ontology annotations summary info")
        go_tab = self.h5.get_node("/Annotations/GeneOntology")
        prot_tab = self.h5.get_node("/Protein/Entries")
        exp_codes = frozenset([b"EXP", b"IDA", b"IPI", b"IMP", b"IGI" b"IEP"])
        cnts = collections.Counter()
        cur_enr = None
        for (enr, term), row_iter in itertools.groupby(
            go_tab, operator.itemgetter("EntryNr", "TermNr")
        ):
            evidences = {row["Evidence"] for row in row_iter}
            is_iea = b"IEA" in evidences
            evidences.discard(b"IEA")
            is_exp = not exp_codes.isdisjoint(evidences)
            is_cur = len(evidences.difference(exp_codes)) > 0
            cnts["annotations_any"] += 1
            if is_exp:
                cnts["annotations_exp"] += 1
            if is_cur:
                cnts["annotations_currated"] += 1
            if is_iea:
                cnts["annotations_iea"] += 1
            if cur_enr != enr:
                e = next(prot_tab.where("EntryNr == {}".format(enr))).fetch_all_fields()
                cnts["proteins_any"] += 1
                if e["OmaGroup"] != 0:
                    cnts["protein_OmaGroup"] += 1
                if len(e["OmaHOG"]) > 0:
                    cnts["protein_HOG"] += 1
                cur_enr = enr
        return cnts

    def count_xref_summary(self):
        self.logger.info("Building cross-ref summary info")
        xref_tab = self.h5.get_node("/XRef")
        prot_tab_iter = iter(self.h5.get_node("/Protein/Entries"))
        source = xref_tab.get_enum("XRefSource")
        trusted = frozenset(
            [
                "UniProtKB/SwissProt",
                "UniProtKB/TrEMBL",
                "RefSeq",
                "EntrezGene",
                "Ensembl Gene",
                "Ensembl Protein",
            ]
        )
        if len(trusted.difference(source._names.keys())) > 0:
            raise ValueError("set of trusted xrefs is invalid")
        cnts = collections.Counter()

        entry = next(prot_tab_iter)
        for enr, xref_it in itertools.groupby(xref_tab, operator.itemgetter("EntryNr")):
            while entry["EntryNr"] < enr:
                entry = next(prot_tab_iter)
            sources_all = [source._values[x["XRefSource"]] for x in xref_it]
            cnts += collections.Counter(sources_all)
            has_trusted_xref = len(trusted.intersection(sources_all)) > 0
            if has_trusted_xref:
                cnts["trusted_all"] += 1
                if entry["OmaGroup"] != 0:
                    cnts["trusted_OmaGroup"] += 1
                if len(entry["OmaHOG"]) > 0:
                    cnts["trusted_HOG"] += 1
        return cnts

    def collect_group_sizes(self):
        self.logger.info("Building grouping size histograms")
        groupings = ("OmaHOG", "OmaGroup")
        memb_cnts = {grp: collections.defaultdict(int) for grp in groupings}
        prot_tab = self.h5.get_node("/Protein/Entries")
        for row in prot_tab:
            for grp in groupings:
                if grp == "OmaHOG":
                    m = hog_re.match(row[grp])
                    if m is None:
                        continue
                    grp_id = int(m.group("fam"))
                else:
                    grp_id = int(row[grp])
                    if grp_id == 0:
                        continue
                memb_cnts[grp][grp_id] += 1
        sizes = {grp: collections.defaultdict(int) for grp in groupings}
        for grp in groupings:
            for grp_size in memb_cnts[grp].values():
                sizes[grp][grp_size] += 1
        return sizes

    def compute_domaincovered_sites(self):
        dom_tab = self.h5.get_node("/Annotations/Domains")
        domains = pandas.DataFrame.from_records(dom_tab[:])

        def dlen(coords):
            doms = [int(pos) for pos in coords.split(b":")]
            return sum((doms[i + 1] - doms[i] + 1 for i in range(0, len(doms), 2)))

        # sum all parts of each domain region and store total length in DLen column
        domains = domains.assign(DLen=domains["Coords"].apply(dlen))
        # sum over all domains per protein
        cov_sites = domains.groupby("EntryNr").agg({"DLen": sum})
        return cov_sites

    def add_domain_covered_sites_counts(self):
        """Stores the number of AA covered by a DomainAnnotation.

        This method adds to the hdf5 file a /Protein/DomainCoverage array that
        contains the number of AA sites covered by a domain. The position
        corresponds to the protein entry numbers in /Protein/Entries.

        :Note: The method assumes that the domains are all non-overlapping.
            If they are not, the reported coverage will be too high!

        :return: covered fractions by domains for each protein
        :rtype: numpy.array"""
        self.logger.info("Counting covered sites by domains")
        cov_sites_df = self.compute_domaincovered_sites()

        prot_tab = self.h5.get_node("/Protein/Entries")
        enr_col = prot_tab.col("EntryNr")
        assert numpy.all(numpy.equal(enr_col, numpy.arange(1, len(prot_tab) + 1)))

        cov_sites = numpy.zeros(len(prot_tab), dtype=numpy.uint32)
        for eNr, coverage in zip(cov_sites_df.index, cov_sites_df.DLen.values):
            cov_sites[eNr - 1] = coverage
        create_node = False
        try:
            dom_cov_tab = self.h5.get_node("/Protein/CoveredSitesByDomains")
            if len(dom_cov_tab) != len(cov_sites):
                self.h5.remove_node("/Protein/CoveredSitesByDomains")
                create_node = True
        except tables.NoSuchNodeError:
            create_node = True
        if create_node:
            dom_cov_tab = self.h5.create_carray(
                "/Protein",
                "CoveredSitesByDomains",
                tables.UInt32Atom(),
                (len(cov_sites),),
            )
        dom_cov_tab[0 : len(cov_sites)] = cov_sites
        return cov_sites / (prot_tab.col("SeqBufferLength") - 1)

    def add_gene_ontology_term_cnts(self):
        """
        Stores the counts of a gene ontology annotation per term in the database.

        This function also includes the implied parental terms.

        :return: hdf5 Table instance with the counts
        """
        cnts, tot_cnts = self.collect_goterm_freqs()
        cnts = sorted(cnts.items())
        tab = self.create_table_if_needed(
            "/Ontologies",
            "GeneOntologyTermCounts",
            description=tablefmt.GeneOntologyTermCounts,
        )
        tab.append(cnts)
        self.h5.set_node_attr(tab, "total_molecular_function", tot_cnts[0])
        self.h5.set_node_attr(tab, "total_biological_process", tot_cnts[1])
        self.h5.set_node_attr(tab, "total_cellular_component", tot_cnts[2])
        return tab

    def add_sequence_suffix_array(self, k=6, fn=None, sa=None):
        """
        Adds the sequence suffix array to the database. NOTE: this
        (obviously) requires A LOT of memory for large DBs.
        """
        # Ensure we're run in correct order...
        assert "Protein" in self.h5.root, "Add proteins before calc. SA!"
        idx_compr = tables.Filters(complevel=6, complib="blosc")

        # Add to separate file if fn is set.
        if fn is None:
            db = self.h5
        else:
            fn = os.path.normpath(
                os.path.join(os.getenv("DARWIN_BROWSERDATA_PATH", ""), fn)
            )
            db = tables.open_file(fn, "w", filters=idx_compr)
            db.create_group("/", "Protein")
            db.root._f_setattr("conversion_start", time.strftime("%c"))
            self.logger.info("opened {}".format(db.filename))

        # Load sequence buffer to memory - this is required to calculate the SA.
        # Do it here (instead of in PySAIS) so that we can use it for computing
        # the split points later.
        seqs = self.h5.get_node("/Protein/SequenceBuffer")[:].tobytes()
        n = len(self.h5.get_node("/Protein/Entries"))

        # Compute & save the suffix array to DB. TODO: work out what compression
        # works best!
        if sa is None:
            sa = sais(seqs)
            sa[:n].sort()  # Sort delimiters by position.
        db.create_carray(
            "/Protein",
            name="SequenceIndex",
            title="concatenated protein sequences suffix array",
            obj=sa,
            filters=idx_compr,
        )

        # Create lookup table for fa2go
        dtype = numpy.uint32 if (n < numpy.iinfo(numpy.uint32).max) else numpy.uint64
        idx = numpy.zeros(sa.shape, dtype=dtype)
        mask = numpy.zeros(sa.shape, dtype=bool)

        # Compute mask and entry index for sequence buff
        for i in range(n):
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
        kmer_lookup_arr = db.create_vlarray(
            "/Protein",
            name="KmerLookup",
            atom=atom(shape=()),
            title="kmer entry lookup table",
            filters=idx_compr,
            expectedrows=len(kmers),
        )
        kmer_lookup_arr._f_setattr("k", k)

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

        if db.filename != self.h5.filename:
            self.logger.info("storing external links to SequenceIndex and KmerLookup")
            self.h5.create_external_link(
                "/Protein",
                "KmerLookup",
                self._relative_path_to_external_node(kmer_lookup_arr),
            )
            self.h5.create_external_link(
                "/Protein",
                "SequenceIndex",
                self._relative_path_to_external_node(db.root.Protein.SequenceIndex),
            )
            db.root._f_setattr("conversion_end", time.strftime("%c"))
            db.close()
            self.logger.info("closed {}".format(db.filename))

    def _relative_path_to_external_node(self, node):
        rel_path = os.path.relpath(
            node._v_file.filename, os.path.dirname(self.h5.filename)
        )
        return str(rel_path + ":" + node._v_pathname)

    def add_hog_domain_prevalence(self):
        # Check that protein entries / domains are added already to the DB
        assert True  # TODO

        # Used later
        hl_tab = self.h5.get_node("/HogLevel")
        create_index_for_columns(hl_tab, "Fam", "IsRoot")
        fam_to_rootlevel = {
            int(row["Fam"]): row["Level"].decode()
            for row in hl_tab.where('~contains(ID, b".") & (IsRoot == True)')
        }

        # Load the HOG -> Entry table to memory
        prot_tab = self.h5.root.Protein.Entries
        # TODO: work out how to do this in a neater way
        df = pandas.DataFrame.from_records(
            (
                (z["EntryNr"], z["OmaHOG"], z["SeqBufferLength"])
                for z in prot_tab.iterrows()
            ),
            columns=["EntryNr", "OmaHOG", "SeqBufferLength"],
        )
        # Strip singletons
        df = df[~(df["OmaHOG"] == b"")]

        # Reformat HOG ID to plain-integer for top-level grouping only
        def root_hog_nr(id_):
            m = hog_re.match(id_)
            return int(m.group("fam"))

        df["OmaHOG"] = df["OmaHOG"].apply(root_hog_nr)

        # Load domains
        domains = pandas.DataFrame.from_records(self.h5.root.Annotations.Domains[:])

        # Ensure sorted by coordinate - TODO: move this to DA import function
        domains["start"] = domains["Coords"].apply(lambda c: int(c.split(b":")[0]))
        domains.sort_values(["EntryNr", "start"], inplace=True)
        domains = domains[["EntryNr", "DomainId"]]

        # Merge domains / entry-hog tables. Keep entries with no domains
        # so that we can count the size of the HOGs.
        df = pandas.merge(df, domains, on="EntryNr", how="left")

        # Gather entry-domain for each HOG.
        hog2dom = []
        hog2info = []
        for hog_id, hdf in tqdm(df.groupby("OmaHOG")):
            size = len(set(hdf["EntryNr"]))

            hdf = hdf[~hdf["DomainId"].isnull()]
            cov = len(set(hdf["EntryNr"]))  # Coverage with any DA

            if (size > 2) and (cov > 1):
                # There are some annotations
                da = collections.defaultdict(list)
                for enum, edf in hdf.groupby("EntryNr"):
                    d = edf["DomainId"]
                    d = tuple(d) if (type(d) != bytes) else (d,)
                    da[d].append(enum)

                da = sorted(da.items(), key=lambda i: len(i[1]), reverse=True)
                c = len(da[0][1])  # Count of prev. DA
                if c > 1:
                    # DA exists in more than one member.
                    cons_da = da[0][0]
                    repr_entry = da[0][1][0]
                    try:
                        tl = hl_tab.read_where("Fam == {}".format(hog_id))[0][
                            "Level"
                        ].decode("ascii")
                    except IndexError as exc:
                        self.logger.exception(
                            "cannot get level for fam {}".format(hog_id)
                        )
                        self.warning(hl_tab.read_where("Fam == {}".format(hog_id)))
                        raise

                    rep_len = hdf[hdf["EntryNr"] == repr_entry]["SeqBufferLength"]
                    rep_len = int(rep_len if len(rep_len) == 1 else list(rep_len)[0])

                    # Save the consensus DA
                    off = len(hog2info)  # Offset in the information table.
                    hog2dom += [(off, d) for d in cons_da]

                    # Save required information about this group for the web
                    # view.
                    hog2info.append(
                        (
                            hog_id,  # HOG ID
                            repr_entry,  # Repr. entry
                            rep_len,  # Repr. entry length
                            tl,  # Top level of HOG
                            size,  # HOG size
                            c,
                        )
                    )  # Prevalence

        # Create tables in file -- done this way as these end up being pretty
        # small tables (<25MB)
        tab = self.h5.create_table(
            "/HOGAnnotations",
            "DomainArchPrevalence",
            tablefmt.HOGDomainArchPrevalenceTable,
            createparents=True,
            expectedrows=len(hog2info),
        )
        self._write_to_table(tab, hog2info)
        tab.flush()  # Required?

        # HOG <-> Domain table
        tab = self.h5.create_table(
            "/HOGAnnotations",
            "Domains",
            tablefmt.HOGDomainPresenceTable,
            createparents=True,
            expectedrows=len(hog2dom),
        )
        self._write_to_table(tab, hog2dom)
        tab.flush()  # Required?

    def add_per_species_aux_groupdata(self):
        self.logger.info("adding per species auxillary group/hog stats")
        omagroup_aux_adder = PerGenomeOMAGroupAuxDataAdder(self.h5)
        omagroup_aux_adder.load_data()
        omagroup_aux_adder.write()
        self.logger.info("  info for OMA Groups written")

        hog_aux_adder = PerGenomeHOGAuxAdder(self.h5)
        hog_aux_adder.load_data()
        hog_aux_adder.write()
        self.logger.info("  info for HOGs written")


def download_url_if_not_present(url, force_copy=False):
    if url.startswith("file://") and not force_copy:
        fname = url[len("file://") :]
        if os.path.exists(fname):
            common.package_logger.info(
                'using file "{}" directly from source without copying.'.format(url)
            )
            return fname
        common.package_logger.warning("file {} does not exist".format(fname))
        return None
    tmpfolder = os.path.join(
        os.getenv("DARWIN_NETWORK_SCRATCH_PATH", "/tmp"), "Browser", "xref"
    )
    basename = url.split("/")[-1]
    fname = os.path.join(tmpfolder, basename)
    if not os.path.exists(tmpfolder):
        os.makedirs(tmpfolder)
    if not os.path.exists(fname):
        common.package_logger.info("downloading {} into {}".format(url, fname))
        try:
            urllib.request.urlretrieve(url, fname)
        except urllib.request.URLError:
            common.package_logger.warn("cannot download {}".format(url))
    return fname


def iter_domains(url):
    DomainTuple = collections.namedtuple("DomainTuple", ("md5", "id", "coords"))

    fname = download_url_if_not_present(url)
    if fname is None:
        return
    with common.auto_open(fname, "rt") as uncompressed:
        dialect = csv.Sniffer().sniff(uncompressed.read(4096))
        uncompressed.seek(0)
        csv_reader = csv.reader(uncompressed, dialect)
        col_md5, col_id, col_coord = (None,) * 3
        coord_fromat_trans = str.maketrans("-,", "::")

        for lineNr, row in enumerate(csv_reader):
            if col_md5 is None:
                # identify which tuples to use.
                if len(row) >= 9:
                    # representative_proteins format. use columns 5-7
                    col_md5, col_id, col_coord = 4, 5, 6
                elif len(row) == 3:
                    # additionally created ones, minimal format
                    col_md5, col_id, col_coord = 0, 1, 2
                else:
                    raise DataImportError(
                        "Unknown Domain Annotation format in {}".format(
                            uncompressed.filename
                        )
                    )
            try:
                dom = DomainTuple(
                    row[col_md5],
                    row[col_id],
                    row[col_coord].translate(coord_fromat_trans),
                )
                if lineNr < 10:
                    # do some sanity checks on the first few lines
                    if re.match(r"[0-9a-f]{32}$", dom.md5) is None:
                        raise DataImportError(
                            "md5 hash of line {:d} has unexpected values: {}".format(
                                lineNr, dom.md5
                            )
                        )
                    if re.match(r"([1-4]\.\d+\.\d+\.\d+|PF\d+)$", dom.id) is None:
                        raise DataImportError(
                            "Domain-ID of line {:d} has unexpected value: {}".format(
                                lineNr, dom.id
                            )
                        )
                    if re.match(r"\d+:\d+", dom.coords) is None:
                        raise DataImportError(
                            "Domain coordinates in line {:d} has unexpected value: {}".format(
                                lineNr, dom.coords
                            )
                        )
                yield dom
            except Exception:
                common.package_logger.exception(
                    "cannot create tuple from line {}".format(lineNr)
                )


def only_pfam_or_cath_domains(iterable):
    cath_re = re.compile(r"[1-4]\.")
    for dom in iterable:
        if dom.id.startswith("PF") or cath_re.match(dom.id) is not None:
            yield dom


def filter_duplicated_domains(iterable):
    """filter duplicated domain annotations that come from different proteins
    with the exact same sequence."""
    seen = set([])
    ignored = 0
    for dom in iterable:
        if not dom in seen:
            seen.add(dom)
            yield dom
        else:
            ignored += 1
    common.package_logger.info(
        "skipped {} duplicated domains. {} distinct domains yielded".format(
            ignored, len(seen)
        )
    )


class RootHOGMetaDataLoader(object):
    """RootHOG Meta data extractor.

    This class provides the means to import the Keywords of the RootHOGs
    into the hdf5 database. The data is stored under in the node defined
    by :attr:`meta_data_path`, which defaults to /RootHOG/MetaData.
    """

    keyword_name = "RootHOG_Keywords.drw"
    expected_keys = ["Keywords"]
    tab_description = tablefmt.RootHOGMetaTable
    meta_data_path = "/RootHOG/MetaData"

    def __init__(self, db):
        self.db = db

    def add_data(self):
        common.package_logger.info("adding {}".format(self.meta_data_path))
        nr_groups = self._get_nr_of_groups()
        has_meta_data = self._check_textfiles_avail()
        if has_meta_data:
            data = self._load_data()
            encoded_data = {}
            for key in self.expected_keys:
                encoded_data[key] = [
                    x.encode("utf-8") if isinstance(x, str) else b"-" for x in data[key]
                ]
                if nr_groups != len(encoded_data[key]):
                    raise DataImportError(
                        "nr of oma groups does not match the number of {}".format(key)
                    )
        else:
            common.package_logger.warning(
                "No fingerprint nor keyword information available"
            )
            encoded_data = {}
            for key in self.expected_keys:
                if key == "Fingerprints":
                    encoded_data[key] = [b"n/a"] * nr_groups
                else:
                    encoded_data[key] = [b""] * nr_groups

        grptab, keybuf = self._create_db_objects(nr_groups)
        self._fill_data_into_db(encoded_data, grptab, keybuf)
        if "NrMembers" in grptab.colnames:
            grptab.modify_column(
                column=self._get_group_member_counts(), colname="NrMembers"
            )
        self._create_indexes(grptab)

    def _create_db_objects(self, nrows):
        key_path = os.path.join(os.path.dirname(self.meta_data_path), "KeywordBuffer")
        try:
            self.db.get_node(self.meta_data_path)
            self.db.remove_node(self.meta_data_path)
            self.db.remove_node(key_path)
        except tables.NoSuchNodeError:
            pass
        root, name = self.meta_data_path.rsplit("/", 1)
        grptab = self.db.create_table(
            root, name, self.tab_description, expectedrows=nrows, createparents=True
        )
        buffer = self.db.create_earray(
            root,
            "KeywordBuffer",
            tables.StringAtom(1),
            (0,),
            "concatenated keywords descriptions",
            expectedrows=500 * nrows,
        )
        return grptab, buffer

    def _fill_data_into_db(self, encoded_data, grp_tab, key_buf):
        row = grp_tab.row
        buf_pos = 0
        keywords = encoded_data["Keywords"]
        for i in range(len(keywords)):
            row["FamNr"] = i + 1
            row["KeywordOffset"] = buf_pos
            row["KeywordLength"] = len(keywords[i])
            row.append()
            key = numpy.ndarray(
                (len(keywords[i]),), buffer=keywords[i], dtype=tables.StringAtom(1)
            )
            key_buf.append(key)
            buf_pos += len(keywords[i])
        grp_tab.flush()
        key_buf.flush()

    def _create_indexes(self, grp_tab):
        create_index_for_columns(grp_tab, "FamNr")

    def _load_data(self):
        return callDarwinExport("GetRootHOGData()")

    def _get_nr_of_groups(self):
        tab = self.db.get_node("/HogLevel")
        try:
            return tab[tab.colindexes["Fam"][-1]]["Fam"]
        except KeyError:
            return max(tab.col("Fam"))

    def _get_group_member_counts(self):
        pass

    def _check_textfiles_avail(self):
        rootdir = os.getenv("DARWIN_BROWSERDATA_PATH", "")
        fn = os.path.join(rootdir, self.keyword_name)
        return os.path.exists(fn)


class OmaGroupMetadataLoader(RootHOGMetaDataLoader):
    """OMA Group Meta data extractor.

    This class provides the means to import the Keywords and Fingerprints
    of the OMA Groups into the hdf5 database. The data is stored under
    in the node defined by :attr:`meta_data_path`, which defaults to
    /OmaGroups/MetaData.
    """

    keyword_name = "Keywords.drw"
    finger_name = "Fingerprints"
    expected_keys = ["Keywords", "Fingerprints"]
    tab_description = tablefmt.OmaGroupTable
    meta_data_path = "/OmaGroups/MetaData"

    def _fill_data_into_db(self, encoded_data, grp_tab, key_buf):
        row = grp_tab.row
        buf_pos = 0
        stable_ids = encoded_data["Fingerprints"]
        keywords = encoded_data["Keywords"]
        for i in range(len(stable_ids)):
            row["GroupNr"] = i + 1
            row["Fingerprint"] = stable_ids[i]
            row["KeywordOffset"] = buf_pos
            row["KeywordLength"] = len(keywords[i])
            row.append()
            key = numpy.ndarray(
                (len(keywords[i]),), buffer=keywords[i], dtype=tables.StringAtom(1)
            )
            key_buf.append(key)
            buf_pos += len(keywords[i])
        grp_tab.flush()
        key_buf.flush()

    def _create_indexes(self, grp_tab):
        create_index_for_columns(grp_tab, "Fingerprint", "GroupNr")

    def _load_data(self):
        return callDarwinExport("GetGroupData()")

    def _get_nr_of_groups(self):
        etab = self.db.get_node("/Protein/Entries")
        try:
            return etab[etab.colindexes["OmaGroup"][-1]]["OmaGroup"]
        except KeyError:
            return max(etab.col("OmaGroup"))

    def _get_group_member_counts(self):
        grp_nr, cnts = numpy.unique(
            self.db.get_node("/Protein/Entries").col("OmaGroup"), return_counts=True
        )
        if grp_nr[0] == 0:
            cnts = cnts[1:]
        assert len(cnts) == self._get_nr_of_groups()
        return cnts

    def _check_textfiles_avail(self):
        rootdir = os.getenv("DARWIN_BROWSERDATA_PATH", "")
        fn1 = os.path.join(rootdir, self.keyword_name)
        fn2 = os.path.join(rootdir, self.finger_name)
        return os.path.exists(fn1) and os.path.exists(fn2)


class DescriptionManager(object):
    def __init__(self, db, entry_path, buffer_path):
        self.db = db
        self.entry_path = entry_path
        self.buffer_path = buffer_path

    def __enter__(self):
        self.entry_tab = self.db.get_node(self.entry_path)
        if not numpy.all(
            numpy.equal(
                self.entry_tab.col("EntryNr"), numpy.arange(1, len(self.entry_tab) + 1)
            )
        ):
            raise RuntimeError("entry table is not sorted")

        root, name = os.path.split(self.buffer_path)
        self.desc_buf = self.db.create_earray(
            root,
            name,
            tables.StringAtom(1),
            (0,),
            "concatenated protein descriptions",
            expectedrows=len(self.entry_tab) * 100,
        )
        self.cur_eNr = None
        self.cur_desc = []
        bufindex_dtype = numpy.dtype(
            [
                (col, self.entry_tab.coldtypes[col])
                for col in ("DescriptionOffset", "DescriptionLength")
            ]
        )
        # columns to be stored in entry table with buffer index data
        self.buf_index = numpy.zeros(len(self.entry_tab), dtype=bufindex_dtype)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.cur_eNr:
            self._store_description()
        self.desc_buf.flush()
        self.entry_tab.modify_columns(
            columns=self.buf_index, names=self.buf_index.dtype.names
        )
        self.entry_tab.flush()

    def add_description(self, eNr, desc):
        """stages a description for addition. Note that the descriptions
        must be ordered according to the entryNr, i.e. all descriptions
        related to eNr X must be staged before changeing to another eNr."""
        if self.cur_eNr and self.cur_eNr != eNr:
            self._store_description()
            self.cur_desc = []
        self.cur_eNr = eNr
        self.cur_desc.append(desc)

    def _store_description(self):
        buf = "; ".join(self.cur_desc).encode("utf-8")
        buf = buf[0 : 2**16 - 1]  # limit to max value of buffer length field
        len_buf = len(buf)
        idx = self.cur_eNr - 1
        self.buf_index[idx]["DescriptionOffset"] = len(self.desc_buf)
        self.buf_index[idx]["DescriptionLength"] = len_buf
        self.desc_buf.append(
            numpy.ndarray((len_buf,), buffer=buf, dtype=tables.StringAtom(1))
        )


class GeneOntologyManager(object):
    ontology_url = "http://purl.obolibrary.org/obo/go/go-basic.obo"

    def __init__(self, db, annotation_path, ontology_path):
        self.db = db
        self.annotation_path = annotation_path
        self.ontology_path = ontology_path
        self._go_buf = []
        self.quote_re = re.compile(r"([[,])([\w_:]+)([,\]])")

    def __enter__(self):
        go_obo_file = download_url_if_not_present(self.ontology_url)
        # check that ontology file is not broken. if we can build it, it should be ok
        self.go = GeneOntology(OntologyParser(go_obo_file))
        self.go.parse()

        with open(go_obo_file, "rb") as fh:
            go_obo = fh.read()
        root, name = os.path.split(self.ontology_path)
        obo = self.db.create_carray(
            root,
            name,
            title="Gene ontology hierarchy definition",
            createparents=True,
            obj=numpy.ndarray(len(go_obo), buffer=go_obo, dtype=tables.StringAtom(1)),
        )
        obo._f_setattr("ontology_release", self._get_obo_version(obo))

        root, name = os.path.split(self.annotation_path)
        self.go_tab = self.db.create_table(
            root,
            name,
            tablefmt.GeneOntologyTable,
            "Gene Ontology annotations",
            expectedrows=1e8,
            createparents=True,
        )
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._flush_buffers()
        self.go_tab.flush()

    def _get_obo_version(self, obo_arr):
        header = obo_arr[0:1000].tobytes()
        rel_info = re.search(rb"data-version:\s*(?P<version>[\w/_ -]+)", header)
        if rel_info is not None:
            rel_info = rel_info.group("version").decode()
        return rel_info

    def _flush_buffers(self):
        common.package_logger.info("flushing go annotations buffers")
        if len(self._go_buf) > 0:
            self.go_tab.append(self._go_buf)
        self._go_buf = []

    def add_annotations(self, enr, gos):
        """parse go annotations and add them to the go buffer"""
        if not (isinstance(enr, int) and isinstance(gos, str)):
            raise ValueError("input data invalid")
        for t in gos.split("; "):
            t = t.strip()
            try:
                term, rem = t.split("@")
            except ValueError as e:
                common.package_logger.warning("cannot parse GO annotation: " + t)
                continue

            try:
                term_nr = self.go.term_by_id(term).id
            except ValueError:
                common.package_logger.warning(
                    "invalid GO term for entry {:d}: {:s} (likely obsolete)".format(
                        enr, term
                    )
                )
                continue
            rem = rem.replace("{", "[")
            rem = rem.replace("}", "]")
            rem = self.quote_re.sub(r'\g<1>"\g<2>"\g<3>', rem)
            for evi, refs in eval(rem):
                for ref in refs:
                    self._go_buf.append((enr, term_nr, evi, ref.encode("utf-8")))
            if len(self._go_buf) > 2e6:
                self._flush_buffers()


class GroupAnnotatorInclGeneRefs(familyanalyzer.GroupAnnotator):
    def _annotateGroupR(self, node, og, idx=0):
        is_geneRef = familyanalyzer.OrthoXMLQuery.is_geneRef_node(node)
        is_og = self.parser.is_ortholog_group(node)
        if is_og or is_geneRef:
            if self.parser.is_paralog_group(
                node.getparent()
            ) or node.getparent().tag == "{{{ns0}}}groups".format(**self.parser.ns):
                # direct children of a paralogGroup. set IsRoot to true in most general TaxRange property tag
                prop_tag = "{{{ns0}}}property".format(**self.parser.ns)
                for child in node[::-1]:
                    if child.tag == prop_tag and child.attrib["name"] == "TaxRange":
                        child.attrib["IsRoot"] = str(True)
                        break

        if is_geneRef:
            node.set("LOFT", og)
        else:
            super()._annotateGroupR(node, og, idx)


class HogConverter(object):
    def __init__(self, entry_tab, release_char=None, tax_tab=None, tax_2_code=None):
        self.fam_re = re.compile(r"HOG:(?P<release>[A-Z]+)?(?P<fam_nr>\d+)")
        self.hogs = numpy.zeros(
            shape=(len(entry_tab) + 1,), dtype=entry_tab.cols.OmaHOG.dtype
        )
        self.entry_tab = entry_tab
        self.taxrange_2_taxid = (
            self._extract_taxrange_2_taxid_map(tax_tab) if tax_tab else None
        )
        self.taxid_2_code = tax_2_code
        if release_char is None:
            self.release_char = ""
        elif re.match(r"^[A-Z]?$", release_char):
            self.release_char = release_char
        else:
            raise ValueError(
                "invalid release_char value: {}. Expected is a single capital ascii character".format(
                    release_char
                )
            )

    def _extract_taxrange_2_taxid_map(self, tax_tab):
        return {row["Name"].decode(): int(row["NCBITaxonId"]) for row in tax_tab}

    def attach_newick_taxonomy(self, tree):
        self.taxonomy = familyanalyzer.NewickTaxonomy(tree)

    def _assert_hogid_has_correct_prefix(self, fa_parser):
        for grp in fa_parser.getToplevelGroups():
            if not grp.get("id").startswith("HOG:{:s}".format(self.release_char)):
                grp.set(
                    "id", "HOG:{:s}{:07d}".format(self.release_char, int(grp.get("id")))
                )

    def convert_file(self, fn, store=None):
        p = familyanalyzer.OrthoXMLParser(fn)
        self._assert_hogid_has_correct_prefix(p)
        if hasattr(self, "taxonomy"):
            tax = self.taxonomy
        else:
            tax = familyanalyzer.TaxonomyFactory.newTaxonomy(p)
        annotator = GroupAnnotatorInclGeneRefs(
            p, self.taxrange_2_taxid, self.taxid_2_code
        )
        annotator.annotateMissingTaxRanges(tax)
        annotator.annotateDoc()

        levs = []
        for fam in p.getToplevelGroups():
            m = self.fam_re.match(self.get_hog_id(fam))
            fam_nr = int(m.group("fam_nr"))
            for taxnode in familyanalyzer.OrthoXMLQuery.getTaxRangeNodes(
                fam, recursively=True
            ):
                ognode = taxnode.getparent()
                levs.append(
                    (fam_nr, self.get_hog_id(ognode), taxnode.get("value"))
                    + self.get_hog_scores(ognode, taxnode)
                    + (
                        self.get_nr_member_genes(ognode),
                        bool(taxnode.get("IsRoot", False)),
                        -1,
                    )
                )

        geneNodes = p.root.findall(
            ".//{{{ns0}}}geneRef".format(**familyanalyzer.OrthoXMLParser.ns)
        )
        for x in geneNodes:
            self.hogs[int(x.get("id"))] = x.get("LOFT")
        if store is not None:
            p.write(store, pretty_print=True)
        return levs

    def write_hogs(self):
        """update the Entry Table with the newly collected OmaHOG values for all
        the proteins at once.

        .. note: This method will overwrite any previous value of the OmaHOG column"""
        self.entry_tab.modify_column(0, len(self.entry_tab), 1, self.hogs[1:], "OmaHOG")
        self.entry_tab.flush()

    def get_hog_id(self, node):
        try:
            id = node.get("id")
            id = id.split("_")[0]
        except AttributeError:
            id = node.get("og")
        return id

    def get_hog_scores(self, og_node, tax_node):
        """extract the scores associated with an orthologGroup node

        only scores that are defined in HOGsTable are extract. The method
        returns a tuple with the scores in the order of the score fields."""
        all_score_ids = ("CompletenessScore", "ImpliedLosses")
        parse_fun = {"CompletenessScore": float, "ImpliedLosses": int}
        scores = collections.OrderedDict(
            [(score, tablefmt.HOGsTable.columns[score].dflt) for score in all_score_ids]
        )
        for score in og_node.iterfind("{*}score"):
            score_id = score.get("id")
            scores[score_id] = parse_fun[score_id](score.get("value"))
        # might be overwritten in tax_node (more specific if available)
        for score_id in all_score_ids:
            val = tax_node.get(score_id)
            if val is not None:
                scores[score_id] = parse_fun[score_id](val)
        return tuple(scores.values())

    def get_nr_member_genes(self, og_node):
        for child in og_node:
            if (
                child.tag
                == "{{{ns0}}}property".format(**familyanalyzer.OrthoXMLParser.ns)
                and child.get("name") == "NrMemberGenes"
            ):
                return int(child.get("value"))
        common.package_logger.warning(
            "couldn't find NrMemberGenes property. scanning xml file for geneRefs"
        )
        return len(familyanalyzer.OrthoXMLQuery.getGeneRefNodes(og_node))


class XRefImporter(object):
    """Object to import various types of crossreferences into hdf5.

    The XRefImporter registers at a db_parser object various handlers
    to import the various types of xrefs, namely ids, go-terms,
    EC annotations and descriptions."""

    def __init__(
        self,
        db_parser,
        genomes_tab,
        xref_tab,
        ec_tab,
        go_manager,
        desc_manager,
        approx_adder,
        uniprot_xrefs_adder,
    ):
        self.xrefs = []
        self.ec = []
        self.xref_tab = xref_tab
        self.ec_tab = ec_tab
        self.go_manager = go_manager
        self.desc_manager = desc_manager
        self.approx_adder = approx_adder
        self.uniprot_xref_adder = uniprot_xrefs_adder

        self.verif_enum = tablefmt.XRefTable.columns.get("Verification").enum
        xrefEnum = tablefmt.XRefTable.columns.get("XRefSource").enum
        tag_to_enums = {
            "GI": (xrefEnum["GI"], "exact"),
            "EntrezGene": (xrefEnum["EntrezGene"], "exact"),
            "WikiGene": (xrefEnum["WikiGene"], "unchecked"),
            "IPI": (xrefEnum["IPI"], "unchecked"),
            "Refseq_ID": (xrefEnum["RefSeq"], "exact"),
            "SwissProt": (xrefEnum["UniProtKB/SwissProt"], "exact"),
            "GeneName": (xrefEnum["Gene Name"], "unchecked"),
            "ORFNames": (xrefEnum["ORF Name"], "unchecked"),
            "OrderedLocusNames": (xrefEnum["Ordered Locus Name"], "unchecked"),
            "ProtName": (xrefEnum["Protein Name"], "unchecked"),
            "Synonyms": (xrefEnum["Synonym"], "unchecked"),
            "HGNC_Id": (xrefEnum["HGNC"], "unchecked"),
            "SMR": (xrefEnum["Swiss Model"], "unchecked"),
            "PDB": (xrefEnum["PDB"], "unchecked"),
            "EMBL": (xrefEnum["EMBL"], "unchecked"),
            "ID": (xrefEnum["SourceID"], "exact"),
            "AC": (xrefEnum["SourceAC"], "exact"),
        }
        for tag, enumval in tag_to_enums.items():
            db_parser.add_tag_handler(
                tag,
                lambda key, enr, typ=enumval: self.multi_key_handler(
                    key, enr, typ[0], typ[1]
                ),
            )
        db_parser.add_tag_handler(
            "DE", lambda key, enr: self.description_handler(key, enr)
        )
        db_parser.add_tag_handler("GO", self.go_handler)
        db_parser.add_tag_handler("ID", self.assign_source_handler)
        db_parser.add_tag_handler("AC", self.assign_source_handler)
        db_parser.add_tag_handler("EC", self.ec_handler)

        for tag in ["SwissProt_AC", "UniProt"]:  # UniProt/TrEMBL tag is cut to UniProt!
            db_parser.add_tag_handler(
                tag,
                lambda key, enr, typ=xrefEnum[
                    "UniProtKB/TrEMBL"
                ]: self.remove_uniprot_code_handler(key, enr, typ),
            )

        db_parser.add_post_entry_handler(self.add_approx_xrefs)
        # register the potential_flush as end_of_entry_notifier
        db_parser.add_post_entry_handler(self.potential_flush)

        self.db_parser = db_parser
        self.xrefEnum = xrefEnum
        self.ENS_RE = re.compile(
            r"ENS(?P<species>[A-Z]{0,3})(?P<typ>[GTP])(?P<num>\d{11})"
        )
        self.FB_RE = re.compile(r"FB(?P<typ>[gnptr]{2})(?P<num>\d{7})")
        self.NCBI_RE = re.compile(r"[A-Z]{3}\d{5}\.\d$")
        self.WB_RE = re.compile(r"WBGene\d{8}$")
        self.EC_RE = re.compile(r"\d+\.(\d+|-)\.(\d+|-)\.(\d+|-)")
        self.ENSGENOME_RE = re.compile(
            b"Ensembl (Metazoa|Plant|Fungi|Protist|Bacteria)", re.IGNORECASE
        )

        self.FLUSH_SIZE = 5e6

        # info about current genome
        self.genomes_tab = genomes_tab
        self._cur_genome = None

    def _get_genome_info(self, entry_nr):
        if not (
            self._cur_genome is not None
            and self._cur_genome["EntryOff"]
            < entry_nr
            <= self._cur_genome["EntryOff"] + self._cur_genome["TotEntries"]
        ):
            self._cur_genome = self.genomes_tab[
                self.genomes_tab["EntryOff"].searchsorted(entry_nr + 1) - 1
            ]
        return self._cur_genome

    def from_EnsemblGenome(self, entry_nr):
        genome_info = self._get_genome_info(entry_nr)
        return self.ENSGENOME_RE.search(genome_info["Release"]) is not None

    def flush_buffers(self):
        common.package_logger.info("flushing xrefs and ec buffers")
        if len(self.xrefs) > 0:
            self.xref_tab.append(sorted(uniq(self.xrefs, transform=lambda x: x[:3])))
            self.xrefs = []
        if len(self.ec) > 0:
            self.ec_tab.append(sorted(uniq(self.ec)))
            self.ec = []

    def potential_flush(self, enr=None):
        if len(self.xrefs) > self.FLUSH_SIZE:
            self.flush_buffers()

    def _add_to_xrefs(self, eNr, enum_nr, key, verif="unchecked"):
        if not isinstance(eNr, int):
            raise ValueError("eNr is of wrong type:" + str(eNr))
        self.xrefs.append((eNr, enum_nr, key.encode("utf-8"), self.verif_enum[verif]))

    def key_value_handler(self, key, eNr, enum_nr, verif="unchecked"):
        """basic handler that simply adds a key (the xref) under a given enum_nr"""
        self._add_to_xrefs(eNr, enum_nr, key, verif)

    def multi_key_handler(self, multikey, eNr, enum_nr, verif="unchecked"):
        """try to split the myltikey field using '; ' as a delimiter and add each
        part individually under the passed enum_nr id type."""
        for key in multikey.split("; "):
            if key.startswith("Rep"):
                continue
            pos = key.find(".Rep")
            if pos > 0:
                key = key[0:pos]
            self._add_to_xrefs(eNr, enum_nr, key, verif)

    def assign_source_handler(self, multikey, eNr):
        """handler that splits the multikey field at '; ' locations and
        tries to guess for each part the id_type. If a type could be
        identified, it is added under with this id type, otherwise left out."""
        for key in multikey.split("; "):
            ens_match = self.ENS_RE.match(key)
            if ens_match is not None:
                typ = ens_match.group("typ")
                if typ == "P":
                    enum_nr = self.xrefEnum["Ensembl Protein"]
                elif typ == "G":
                    enum_nr = self.xrefEnum["Ensembl Gene"]
                elif typ == "T":
                    enum_nr = self.xrefEnum["Ensembl Transcript"]
                common.package_logger.debug(
                    "ensembl: ({}, {}, {})".format(key, typ, enum_nr)
                )
                self._add_to_xrefs(eNr, enum_nr, key, "exact")

            for enum, regex in {
                "FlyBase": self.FB_RE,
                "NCBI": self.NCBI_RE,
                "WormBase": self.WB_RE,
            }.items():
                match = regex.match(key)
                if match is not None:
                    enum_nr = self.xrefEnum[enum]
                    self._add_to_xrefs(eNr, enum_nr, key, "unchecked")
            if self.from_EnsemblGenome(eNr):
                self._add_to_xrefs(eNr, self.xrefEnum.EnsemblGenomes, key, "exact")

    def go_handler(self, gos, enr):
        self.go_manager.add_annotations(enr, gos)

    def ec_handler(self, ecs, enr):
        for t in ecs.split("; "):
            t = t.strip()
            acc_match = self.EC_RE.match(t)
            if acc_match is not None:
                self.ec.append((enr, acc_match.group(0)))

    def description_handler(self, de, eNr):
        self.desc_manager.add_description(eNr, de)

    def remove_uniprot_code_handler(self, multikey, eNr, enum_nr):
        """remove the species part (sep by '_') of a uniprot long accession to the short acc"""
        common.package_logger.debug(
            "remove_uniprot_code_handler called ({}, {},{})".format(
                multikey, eNr, enum_nr
            )
        )
        for key in multikey.split("; "):
            pos = key.find("_")
            store_key = key if pos < 0 else key[:pos]
            self._add_to_xrefs(eNr, enum_nr, store_key, "exact")
            for tup in self.uniprot_xref_adder.iter_xreftuples_for_up(store_key, eNr):
                self.xrefs.append(tup)

    def add_approx_xrefs(self, enr):
        self.xrefs.extend(self.approx_adder.iter_approx_xrefs_for(enr))


class UniProtAdditionalXRefImporter(object):
    def __init__(self, *args):
        self.verify_enum = tablefmt.XRefTable.columns.get("Verification").enum
        xref_enum = tablefmt.XRefTable.columns.get("XRefSource").enum
        self.mapping = {
            x: xref_enum[x]
            for x in (
                "STRING",
                "DisGeNET",
                "neXtProt",
                "Bgee",
                "HGNC",
                "EPD",
                "GlyConnect",
                "ChEMBL",
                "SwissPalm",
            )
        }
        self.mapping.update(
            {
                "GeneID": xref_enum["EntrezGene"],
                "GeneWiki": xref_enum["WikiGene"],
                "SMR": xref_enum["Swiss Model"],
                "Ensembl": xref_enum["Ensembl Transcript"],
            }
        )
        self._lookup = collections.defaultdict(list)
        for fpath in args:
            if os.path.exists(fpath):
                self._load_projected_xrefs(fpath)
            else:
                common.package_logger.warning(
                    'UniProt projected XRef file "{}" does not exist. Skipping'.format(
                        fpath
                    )
                )

    def _load_projected_xrefs(self, fpath):
        with open(fpath, "rt", newline="", encoding="utf-8", errors="ignore") as fh:
            csv_reader = csv.reader(fh, delimiter="\t")
            for row in csv_reader:
                try:
                    source = self.mapping[row[1]]
                    self._lookup[row[0]].append((source, row[2]))
                except KeyError:
                    pass

    def iter_xreftuples_for_up(self, accession, enr, confidence=None):
        if confidence is None:
            confidence = self.verify_enum["unchecked"]
        for xref in self._lookup[accession]:
            yield enr, xref[0], xref[1].encode("utf-8"), confidence


class ApproximateXRefImporter(object):
    def __init__(self, *args):
        self.verif_enum = tablefmt.XRefTable.columns.get("Verification").enum
        xref_enum = tablefmt.XRefTable.columns.get("XRefSource").enum
        self.mapping = {
            "UniProtKB/SwissProt": xref_enum["UniProtKB/SwissProt"],
            "UniProtKB/TrEMBL": xref_enum["UniProtKB/TrEMBL"],
            "EntrezGene": xref_enum["EntrezGene"],
            "Refseq_ID": xref_enum["RefSeq"],
            "OrderedLocusNames": xref_enum["Ordered Locus Name"],
            "ORFNames": xref_enum["ORF Name"],
            "ProtName": xref_enum["Protein Name"],
            "GeneName": xref_enum["Gene Name"],
        }
        self.xrefs = []
        for fpath in args:
            if os.path.exists(fpath):
                self.load_approx_xref_file(fpath)
            else:
                common.package_logger.warning(
                    'Approximate XRef file "{}" does not exist. Skipping'.format(fpath)
                )
        self._pos = 0
        self._last_enr = 0

    def load_approx_xref_file(self, fpath):
        with open(fpath, "rt", newline="", encoding="utf-8", errors="ignore") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if row[1] in self.mapping:
                    enr = int(row[0])
                    self.xrefs.append(
                        (
                            enr,
                            self.mapping[row[1]],
                            row[2].encode("utf-8"),
                            self.verif_enum[row[3]],
                        )
                    )
        self.xrefs.sort(key=lambda x: x[0])
        self._pos, self._last_enr = 0, 0

    def iter_approx_xrefs_for(self, entry_nr):
        if entry_nr < self._last_enr:
            raise InvarianceException(
                "entry_nr should always increase between calls of this method: {} expected > {}".format(
                    entry_nr, self._last_enr
                )
            )
        while self._pos < len(self.xrefs) and self.xrefs[self._pos][0] == entry_nr:
            yield self.xrefs[self._pos]
            self._pos += 1
        self._last_enr = entry_nr


class DarwinDbEntryParser:
    def __init__(self):
        """Initializes a Parser for SGML formatted darwin database file"""
        self.tag_handlers = collections.defaultdict(list)
        self.end_of_entry_notifier = []

    def add_tag_handler(self, tag, handler):
        """add a callback handler for a certain tag"""
        self.tag_handlers[tag].append(handler)
        common.package_logger.debug(
            "# handlers for {}: {}".format(tag, len(self.tag_handlers[tag]))
        )

    def add_post_entry_handler(self, handler):
        self.end_of_entry_notifier.append(handler)

    def parse_entrytags(self, fh):
        """AC, CHR, DE, E, EMBL, EntrezGene, GI, GO, HGNC_Name, HGNC_Sym,
        ID, InterPro, LOC, NR , OG, OS, PMP, Refseq_AC, Refseq_ID, SEQ,
        SwissProt, SwissProt_AC, UniProt/TrEMBL, WikiGene, flybase_transcript_id

        :param fh: an already opened file handle to the darwin database
                   file to be parsed."""
        eNr = 0
        for line in fh:
            line = line.strip()
            if not line.startswith("<E>"):
                common.package_logger.debug("skipping line:" + line)
                continue

            eNr += 1
            common.package_logger.debug(
                "entry {}: {}".format(eNr, line.encode("utf-8"))
            )
            entry = lxml.html.fragment_fromstring(line)
            for tag, handlers in self.tag_handlers.items():
                common.package_logger.debug(
                    "tag {} ({} handlers)".format(tag, len(handlers))
                )
                tag_text = [t.text for t in entry.findall("./" + tag.lower())]
                for value in tag_text:
                    # common.package_logger.debug('value of tag: {}'.format(value.encode('utf-8')))
                    if value is None:
                        continue
                    for handler in handlers:
                        handler(value, eNr)
                        # common.package_logger.debug('called handler {} with ({},{})'.format(
                        #    handler, value.encode('utf-8'), eNr))
            for notifier in self.end_of_entry_notifier:
                notifier(eNr)


DomainDescription = collections.namedtuple(
    "DomainDescription", tables.dtype_from_descr(tablefmt.DomainDescriptionTable).names
)


class CathDomainNameParser(object):
    re_pattern = re.compile(r"(?P<id>[0-9.]*)\s{3,}\w{7}\s{3,}:\s*(?P<desc>.*)")
    source = b"CATH/Gene3D"

    def __init__(self, url):
        self.fname = download_url_if_not_present(url)

    def parse(self):
        with common.auto_open(self.fname, "rt") as fh:
            for line in fh:
                match = self.re_pattern.match(line)
                if match is not None:
                    yield DomainDescription(
                        DomainId=match.group("id").encode("utf-8"),
                        Source=self.source,
                        Description=match.group("desc").encode("utf-8"),
                    )


class PfamDomainNameParser(CathDomainNameParser):
    re_pattern = re.compile(r"(?P<id>\w*)\t\w*\t\w*\t\w*\t(?P<desc>.*)")
    source = b"Pfam"


class PerGenomeOMAGroupAuxDataAdder(object):
    node_postfix = "omagroup"

    def __init__(self, h5: tables.File):
        self.h5 = h5
        self.shared_ogs = None
        self.in_groups = None

    def load_data(self):
        gs = self.h5.get_node("/Genome").read()
        Goff = gs["EntryOff"]
        etab = self.h5.get_node("/Protein/Entries")
        shared_ogs = numpy.zeros((len(gs), len(gs)), dtype="i4")
        in_groups = numpy.zeros(len(gs), dtype="i4")

        def enr_to_gnr(entry_nr):
            return Goff.searchsorted(entry_nr, side="left") - 1

        grp_members = self.get_group_members(etab)
        for grp, memb in tqdm(grp_members.items()):
            t0 = time.time()
            gnrs = collections.Counter(map(enr_to_gnr, memb))
            for gn, cnts in gnrs.items():
                in_groups[gn] += cnts
            for g1, g2 in itertools.combinations(gnrs, 2):
                shared_ogs[g1, g2] += 1
                shared_ogs[g2, g1] += 1
            common.package_logger.info(
                "grp {} with {} members took {:.1f}msec".format(
                    grp, len(memb), 1000 * (time.time() - t0)
                )
            )
        self.shared_ogs = shared_ogs
        self.in_groups = in_groups

    def write(self, node=None):
        root = "/Summary" if node is None else node
        for label, data in (
            ("shared_", self.shared_ogs),
            ("prots_in_", self.in_groups),
        ):
            node_full_name = label + self.node_postfix
            try:
                self.h5.remove_node(root, node_full_name)
            except tables.NoSuchNodeError:
                pass
            self.h5.create_array(root, node_full_name, obj=data)

    def get_group_members(self, etab):
        grp_members = collections.defaultdict(list)
        for row in etab.where("OmaGroup > 0"):
            grp_members[row["OmaGroup"]].append(row["EntryNr"])
        return grp_members


class PerGenomeHOGAuxAdder(PerGenomeOMAGroupAuxDataAdder):
    node_postfix = "hog"

    def get_group_members(self, etab):
        grp_members = collections.defaultdict(list)
        for row in etab.where('OmaHOG != b""'):
            fam = int(hog_re.match(row["OmaHOG"]).group("fam"))
            grp_members[fam].append(row["EntryNr"])
        return grp_members


class InvarianceException(Exception):
    pass


def augment_genomes_json_download_file(fpath, h5, backup=".bak"):
    """Augment the genomes.json file in the download section with additional info

    This function stores the ncbi taxonomy identifiers of internal nodes and adds
    the number of ancestral genes to the internal nodes.

    :param fpath: path to genomes.json file
    :param h5: hdf5 database handle."""
    common.package_logger.info("Augmenting genomes.json file with Nr of HOGs per level")

    # load taxonomy and sorter by Name
    tax = h5.get_node("/Taxonomy").read()
    sorter = numpy.argsort(tax["Name"])
    gs = h5.get_node("/Genome").read()
    pe = h5.get_node("/Protein/Entries")
    with open(fpath, "rt") as fh:
        genomes = json.load(fh)
    os.rename(fpath, fpath + ".bak")

    def traverse(node, parent_hogs=None, parent_hogs_support=None):
        if parent_hogs is None:
            parent_hogs = numpy.array([], dtype=h5.get_node("/HogLevel").dtype)
        if parent_hogs_support is None:
            parent_hogs_support = numpy.array([], dtype=h5.get_node("/HogLevel").dtype)

        try:
            n = node["name"].encode("utf-8")
            idx = numpy.searchsorted(tax["Name"], n, sorter=sorter)
            if tax["Name"][sorter[idx]] == n:
                node["taxid"] = taxid = int(tax["NCBITaxonId"][sorter[idx]])
            elif n == b"LUCA":
                taxid = 0
            elif b"(disambiguate" in n:
                # this is a special case to deal with internal species
                # hoginfo is stored at tax[UniProtSpeciesCode]
                taxid = node["id"]
            else:
                node["nr_hogs"] = 0
                raise ValueError("not in taxonomy: {}".format(n))

            for modif, completeness, parent in zip(
                ["", "_support"], [0.0, 0.2], [parent_hogs, parent_hogs_support]
            ):
                hogs = h5.get_node(f"/AncestralGenomes/tax{taxid}/Hogs").read_where(
                    "CompletenessScore > completeness"
                )
                if modif == "":
                    hog_level = hogs
                else:
                    hog_level_support = hogs
                node["nr_hogs{}".format(modif)] = len(hogs)
                diff_parent, dupl_events = hoghelper.compare_levels(
                    parent, hogs, return_duplication_events=True
                )
                changes = collections.defaultdict(int)
                for x in diff_parent["Event"]:
                    changes[x.decode()] += 1
                changes["duplications"] = dupl_events
                if "children" not in node and modif != "_support":
                    # dealing with an extend species, special way to assess gains, based on
                    # HOG singletons that are main isoforms
                    assert changes["gained"] == 0
                    g = numpy.extract(
                        gs["UniProtSpeciesCode"] == node["id"].encode("utf-8"), gs
                    )[0]
                    query_main_iso_of_genome = "(EntryNr > {}) & (EntryNr <= {}) & ((AltSpliceVariant == 0) | (AltSpliceVariant == EntryNr))".format(
                        g["EntryOff"], g["EntryOff"] + g["TotEntries"]
                    )
                    nr_genes, nr_gains = 0, 0
                    for p in pe.where(query_main_iso_of_genome):
                        nr_genes += 1
                        if p["OmaHOG"] == b"":
                            nr_gains += 1
                    changes["gained"] = nr_gains
                    node["nr_genes{}".format(modif)] = nr_genes
                node["evolutionaryEvents{}".format(modif)] = changes
            common.package_logger.info(
                "node: {}: events: {}; events_support: {}".format(
                    n, node["evolutionaryEvents"], node["evolutionaryEvents_support"]
                )
            )
        except Exception:
            common.package_logger.exception("Cannot identify taxonomy id")
            hog_level = parent_hogs.copy()
            hog_level_support = parent_hogs_support.copy()

        if "children" in node:
            for child in node["children"]:
                traverse(
                    child, parent_hogs=hog_level, parent_hogs_support=hog_level_support
                )

    traverse(genomes)
    with open(fpath, "wt") as fh:
        json.dump(genomes, fh)


def getLogger(level="DEBUG"):
    import logging

    log = logging.getLogger("pyoma")
    if isinstance(level, str):
        level = logging.getLevelName(level.upper())
        if not isinstance(level, int):
            level = logging.DEBUG
    log.setLevel(level)
    logHandler = logging.StreamHandler()
    logHandler.setLevel(level)
    logHandler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    log.addHandler(logHandler)
    return log


def main(
    name="OmaServer.h5",
    k=6,
    idx_name=None,
    domains=None,
    log_level="INFO",
    release=None,
    phases=None,
):
    idx_name = (name + ".idx") if idx_name is None else idx_name
    if phases is None:
        phases = list(range(1, 10))
    if (
        not all((isinstance(z, int) for z in phases))
        or min(phases) < 1
        or max(phases) > 9
    ):
        raise ValueError("phases argument must be integers between 1 and 9")
    phases = uniq(sorted(phases))

    def phase1():
        x.add_version(release_char=release)
        x.add_species_data()
        x.add_orthologs()
        x.add_same_species_relations()
        x.add_proteins()

    def phase2():
        x.add_sequence_suffix_array(k=k, fn=idx_name)

    def phase3():
        x.add_hogs()
        x.add_xrefs()

    def phase4():
        x.add_synteny_scores()
        x.add_homoeology_confidence()

    def phase5():
        if domains is None:
            domainFiles = ["file:///dev/null"]
        else:
            domainFiles = list(
                map(
                    lambda url: "file://" + url if url.startswith("/") else url, domains
                )
            )
        log.info("loading domain annotations from {}".format(domainFiles))
        x.add_domain_info(
            filter_duplicated_domains(
                only_pfam_or_cath_domains(
                    itertools.chain.from_iterable(map(iter_domains, domainFiles))
                )
            )
        )
        x.add_domainname_info(
            itertools.chain(
                CathDomainNameParser(
                    "http://download.cathdb.info/cath/releases/latest-release/"
                    "cath-classification-data/cath-names.txt"
                ).parse(),
                PfamDomainNameParser(
                    "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
                ).parse(),
            )
        )

    def phase6():
        x.add_canonical_id()
        x.add_group_metadata()
        x.add_hog_domain_prevalence()
        x.add_roothog_metadata()
        x.add_gene_ontology_term_cnts()

    def phase7():
        x.create_indexes()
        x.update_summary_stats()
        x.add_per_species_aux_groupdata()

    def phase8():
        x.add_cache_of_hogs_by_level(min(os.cpu_count(), 32))

    def phase9():
        genomes_json_fname = os.path.normpath(
            os.path.join(
                os.path.dirname(x.h5.filename), "..", "downloads", "genomes.json"
            )
        )
        augment_genomes_json_download_file(genomes_json_fname, x.h5)

    log = getLogger(log_level)
    log.info("Running the following phases: {}".format(phases))
    x = DarwinExporter(name, logger=log)
    for phase_nr in phases:
        phase_func = locals()["phase{}".format(phase_nr)]
        log.info("calling phase{}".format(phase_nr))
        phase_func()
        log.info("done with phase{}".format(phase_nr))
    x.close()
