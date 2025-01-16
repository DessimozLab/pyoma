import collections
import itertools
import json
import logging
import pickle
import sys
import warnings
from os.path import exists, getsize, join
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from shutil import copy

import ete3
import pandas
import tables
from tables import PerformanceWarning
import omataxonomy


from .builder import DBBuilder, OmaGroupsProvider, XrefStorer
from . import hogconvert
from . import xref as xref_build
from .cachebuilder import create_job_files, combine_results, process_job_file
from ..convert import (
    iter_domains,
    filter_duplicated_domains,
    only_pfam_or_cath_domains,
    CathDomainNameParser,
    PfamDomainNameParser,
    augment_genomes_json_download_file,
)
from ..models import Genome
from ..xref_contrib import reduce_xrefs
from ...common import auto_open

logger = logging.getLogger(__name__)


def phase_genomes(conf):
    with DBBuilder(conf.db, mode="write", logger=logger) as db:
        db.add_version(conf.rel_char)
        db.add_taxonomy(conf.tax_tsv)
        db.add_species_data(conf.gs_tsv)
        with XrefStorer(conf.xref_db, index_cols=["EntryNr"]) as xref_storer:
            db.add_proteins(conf.genomes, OmaGroupsProvider(conf.oma_groups), xref_collector=xref_storer)


def phase_build_seq_indexes(conf):
    with DBBuilder(conf.db, mode="read", logger=logger) as db, DBBuilder(
        conf.out, mode="write", logger=logger, complib="blosc"
    ) as out:
        seqs = db.h5.get_node("/Protein/SequenceBuffer").read().tobytes()
        nr_entries = len(db.h5.get_node("/Protein/Entries"))
        out.add_sequence_index(seqs=seqs, nr_entries=nr_entries, k=6)


def phase_convert_hogs(conf):
    parser = hogconvert.get_orthoxml_parser(conf.db, conf.oma_prot_id)
    if conf.orthoxml_out is not None:
        hogconvert.FullOrthoXMLObserver(parser, conf.orthoxml_out)
    if conf.augmented_orthoxml_out is not None:
        hogconvert.FullAugmentedOrthoXMLObserver(parser, conf.augmented_orthoxml_out)
    with hogconvert.HOGtoHDF5(parser, conf.hdf5_out) as hog_h5:
        annotator = hogconvert.Annotator(parser)
        hogconvert.PerFamilyHOGObserver(parser, hog_h5.store_orthoxml, hog_h5.store_orthoxml_augmented)
        hogconvert.parse_orthoxml(conf.orthoxml, parser)
    # add the per-level HOG tables to the HDF5
    with DBBuilder(conf.hdf5_out, mode="append", logger=logger) as builder:
        with tables.open_file(conf.db) as taxdb:
            lev2tax = {row["Name"]: int(row["NCBITaxonId"]) for row in taxdb.get_node("/Taxonomy").read()}
        builder.add_cache_of_hogs_by_level(lev2tax=lev2tax, nr_procs=4)


def phase_vps(conf):
    with tables.open_file(conf.db, "r") as db:
        with DBBuilder(conf.hdf5_out, mode="write", logger=logger, complib="blosc") as out:
            genomes_tab = db.get_node("/Genome")
            out.add_orthologs(conf.vps_base, genomes=genomes_tab)


def phase_add_domains(conf):
    with tables.open_file(conf.db, "r") as db:
        with DBBuilder(conf.hdf5_out, mode="write", logger=logger) as out:
            md5_to_enr = collections.defaultdict(list)
            for e in db.get_node("/Protein/Entries"):
                md5_to_enr[e["MD5ProteinHash"]].append(e["EntryNr"])
            logger.info("loaded mapping of md5 hashes to entries with %d unique hashes", len(md5_to_enr))
            out.add_domain_info(
                filter_duplicated_domains(
                    only_pfam_or_cath_domains(itertools.chain.from_iterable(map(iter_domains, conf.domains)))
                ),
                md5_to_enr=md5_to_enr,
            )
            out.add_domainname_info(
                itertools.chain(
                    CathDomainNameParser(conf.cath_names).parse(),
                    PfamDomainNameParser(conf.pfam_names).parse(),
                )
            )


def phase_select_alt_splice_variants(conf):
    with DBBuilder(conf.db, mode="append", logger=logger) as db:
        prot_hogid_arr = db.h5.get_node("/OmaHOG")
        hogids = prot_hogid_arr.read()
        db.add_protein_hog_ids(hogids)
        prot_hogid_arr.remove()
        db.identify_and_store_splice_variants(conf.splice_json)


def cache_build_job_generator(conf):
    create_job_files(conf.db, conf.out_prefix)


def cache_build_process_job(conf):
    process_job_file(conf.job_file, conf.db, conf.out)


def cache_build_combine(conf):
    combine_results(job_results=conf.jobs, out=conf.out)


def fetch_refseq(conf):
    xref_build.fetch(
        host=conf.host,
        directory=conf.directory,
        pattern=conf.pattern,
        checksum_ftp_path=conf.checksum_ftp_path,
        out_dir=conf.out,
        ftp_config=conf.ftp_config,
        nr_cpu=conf.nr_cpu,
    )


def build_relevant_taxid_mapping(conf):
    gs = pandas.read_csv(conf.gs_tsv, sep="\t")
    ncbi_taxids = set(gs["OriginalNCBITaxonId"])
    relevant_taxids = xref_build.load_relevant_taxids(ncbi_taxids, omataxonomy.Taxonomy(conf.tax_sqlite))
    with tables.open_file(conf.db, "r") as db:
        genome = pandas.DataFrame(db.get_node("/Genome").read())
        genome["UniProtSpeciesCode"] = genome["UniProtSpeciesCode"].apply(bytes.decode)
    gs = gs.set_index("UniProtSpeciesCode").join(genome.set_index("UniProtSpeciesCode"), how="inner", rsuffix="_db")
    ncbi2oma = collections.defaultdict(set)
    for ncbi, oma in zip(gs["OriginalNCBITaxonId"], gs["NCBITaxonId_db"]):
        ncbi2oma[ncbi].add(oma)

    mapped = {}
    for ncbi, rel_genome_ncbi_tax in relevant_taxids.items():
        genome_ncbis = set()
        for k in rel_genome_ncbi_tax:
            genome_ncbis |= ncbi2oma[k]
        mapped[ncbi] = genome_ncbis
    with auto_open(conf.out, "wb") as fh:
        pickle.dump(mapped, fh)


def filter_and_split_xrefs(conf):
    with auto_open(conf.tax_map, "rb") as fh:
        relevant_taxid_map = pickle.load(fh)
    relevant_taxids = set(relevant_taxid_map.keys())
    with xref_build.ChunkWriter(conf.out_prefix, 30_000) as writer:
        for xref_file in conf.xref:
            xref_build.filter_records_on_taxids(xref_file, writer, relevant_taxids, conf.format)


def map_xrefs(conf):
    with auto_open(conf.tax_map, "rb") as fh:
        relevant_taxid_map = pickle.load(fh)
    xref_build.map_xrefs(
        fpath=conf.xref,
        format=conf.format,
        source=conf.source,
        out_fpath=conf.out,
        db=conf.db,
        seq_idx=conf.seq_idx_db,
        xref_db=conf.xref_source_db,
        taxid_mapping=relevant_taxid_map,
    )


def collect_xrefs(conf):
    xref_build.collect_crossrefs(
        xrefs=conf.xrefs, map_files=conf.map_results, format=conf.format, source=conf.source, out=conf.out
    )


def combine_xrefs(conf):
    xref_build.combine_xrefs(xrefs=conf.xrefs, out=conf.out)


def build_reduced_xrefs(conf):
    reduce_xrefs(conf.db, conf.xrefs, conf.out, nr_procs=conf.nr_procs)


def import_go(conf):
    with auto_open(conf.tax_map, "rb") as fh:
        relevant_taxid_map = pickle.load(fh)
    rel_taxids = set(relevant_taxid_map.keys())
    xref_build.import_go(obo=conf.obo, gafs=conf.gaf, xref_db=conf.xref_db, relevant_taxid=rel_taxids, out=conf.out)


def store_summary_info(conf):
    with DBBuilder(conf.db) as db:
        db.update_summary_stats()
        db.add_hog_domain_prevalence()
        db.add_group_metadata()
        db.add_roothog_metadata()


def gen_aux_files(conf):
    from ..db import Database

    def genome_2_flatgenome_json_dict(g: Genome):
        return {
            "id": g.uniprot_species_code,
            "name": g.sciname,
            "last_updated": g.modification_date("%b %d, %Y"),
            "nr_proteins": g.nr_entries,
            "source": g.release,
            "taxid": g.ncbi_taxon_id,
        }

    with Database(conf.db) as db:
        with auto_open(join(conf.out_dir, "genomes.json"), "wt") as fh:
            tax_without_skiplevel = db.tax.get_subtaxonomy_rooted_at(db.tax.root["NCBITaxonId"], collapse=True)
            genomes = tax_without_skiplevel.as_dict()
            json.dump(genomes, fh)
        augment_genomes_json_download_file(join(conf.out_dir, "genomes.json"), db.get_hdf5_handle())
        with auto_open(join(conf.out_dir, "speciestree.nwk"), "wt") as fh:
            fh.write(db.tax.newick(leaf="sciname", internal="sciname", quoted=True))
        with auto_open(join(conf.out_dir, "speciestree.phyloxml"), "wb") as fh:
            fh.write(db.tax.as_phyloxml())
        with auto_open(join(conf.out_dir, "flatgenomes.json"), "wt") as fh:
            genomes_list = [genome_2_flatgenome_json_dict(g) for g in db.tax.genomes.values()]
            json.dump(genomes_list, fh)


def parse_command_line_args():
    parser = ArgumentParser(description="Builder for OMA Browser hdf5")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity")

    subparsers = parser.add_subparsers(title="Commands")

    genomes_parser = subparsers.add_parser(
        "genomes",
        help="Adding genomes and protein data",
        description="genomes",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    genomes_parser.set_defaults(func=phase_genomes)
    genomes_parser.add_argument("--db", required=True, help="Path to database")
    genomes_parser.add_argument("--gs-tsv", required=True, help="Path to genomes summary file in TSV format")
    genomes_parser.add_argument("--tax-tsv", required=True, help="Path to taxonomy file in TSV format")
    genomes_parser.add_argument("--oma-groups", required=False, help="Path to OMA groups json file")
    genomes_parser.add_argument("--rel-char", required=False, default=None, help="Release character")
    genomes_parser.add_argument("--release", required=False, help="Release of database")
    genomes_parser.add_argument("--xref-db", required=False, help="Path where source xrefs are stored")
    genomes_parser.add_argument(
        "--genomes", required=True, nargs="+", help="List of genome files (json) containing essential data"
    )

    seqindex_parser = subparsers.add_parser("seqindex", help="Adding sequence indexes")
    seqindex_parser.set_defaults(func=phase_build_seq_indexes)
    seqindex_parser.add_argument("--db", required=True, help="Path to database containing sequence")
    seqindex_parser.add_argument("--out", required=True, help="Path to output sequence index database file")

    hogconv_parser = subparsers.add_parser(
        name="hog", help="Converting input orthoxml into reformatted versions and HDF5"
    )
    hogconv_parser.set_defaults(func=phase_convert_hogs)
    hogconv_parser.add_argument("--orthoxml", required=True, help="Path to input orthoxml file")
    hogconv_parser.add_argument("--db", required=True, help="Path to hdf5 database with entires and taxonomy")
    hogconv_parser.add_argument(
        "--hdf5-out", required=True, help="Path to store hoglevel table and per family orthoxml"
    )
    hogconv_parser.add_argument(
        "--augmented-orthoxml-out",
        required=False,
        help="Path where to store the augmented orthoxml file for all the HOGs",
    )
    hogconv_parser.add_argument(
        "--orthoxml-out",
        required=False,
        help="Path where to store the orthoxml file for all the HOGs with updated IDs etc",
    )
    hogconv_parser.add_argument(
        "--oma-prot-id", action="store_true", help="Whether the protId attributes in the input orthoxml contain OMA-IDs"
    )

    vp_parser = subparsers.add_parser("vps", help="Adding pairwise orthologs")
    vp_parser.set_defaults(func=phase_vps)
    vp_parser.add_argument("--db", required=True, help="Path to hdf5 database containing genomes")
    vp_parser.add_argument("--vps-base", required=False, help="Folder where all the pairwise orthologs are stored")
    vp_parser.add_argument("--hdf5-out", required=True, help="Path to store pairwise orthologs in HDF5")

    domain_parser = subparsers.add_parser(
        "domains", help="Adding domain annotations to OMA", formatter_class=ArgumentDefaultsHelpFormatter
    )
    domain_parser.set_defaults(func=phase_add_domains)
    domain_parser.add_argument("--db", required=True, help="Path to database")
    domain_parser.add_argument("--hdf5-out", required=True, help="Path to store domain annotations in HDF5")
    domain_parser.add_argument(
        "--domains",
        required=True,
        nargs="+",
        help="List filenames containing domain annotations for protein sequence hashes",
    )
    domain_parser.add_argument(
        "--cath-names",
        required=False,
        default="http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-names.txt",
        help="Path pointing to cath_names.txt files which provide names for cath domain numbers",
    )
    domain_parser.add_argument(
        "--pfam-names",
        required=False,
        default="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz",
        help="Path pointing to pfam domain name mapping file",
    )

    splice_parser = subparsers.add_parser("splice", help="Adding alternative splice information - set main variant")
    splice_parser.set_defaults(func=phase_select_alt_splice_variants)
    splice_parser.add_argument("--db", required=True, help="Path to database - will be modified")
    splice_parser.add_argument("--splice-json", required=True, help="Path to splice json file")

    # cache builder commands
    cache_job_parser = subparsers.add_parser("cache-job", help="Generate job files for cache building")
    cache_job_parser.set_defaults(func=cache_build_job_generator)
    cache_job_parser.add_argument("--db", required=True, help="Path to database")
    cache_job_parser.add_argument("--out-prefix", default="cache-job", help="prefix for the job files path")

    cache_build_parser = subparsers.add_parser("cache-build", help="Runs cache building for a job-file")
    cache_build_parser.set_defaults(func=cache_build_process_job)
    cache_build_parser.add_argument("--db", required=True, help="Path to database")
    cache_build_parser.add_argument("--job-file", required=True, help="job file to be processed")
    cache_build_parser.add_argument("--out", required=True, help="output filename")

    cache_combine_parser = subparsers.add_parser("cache-combine", help="Combines all cache jobs in a single hdf5")
    cache_combine_parser.set_defaults(func=cache_build_combine)
    cache_combine_parser.add_argument("--jobs", nargs="+", help="job result files to combine")
    cache_combine_parser.add_argument("--out", required=True, help="output path of combined results")

    # refseq subparser
    refseq_fetch_parser = subparsers.add_parser(
        "fetch-refseq", help="Fetching refseq data from remote servers to be integrated"
    )
    refseq_fetch_parser.set_defaults(func=fetch_refseq)
    refseq_fetch_parser.add_argument(
        "--host", default="ftp.ncbi.nih.gov", help="ftp host from where to fetch_parser the data"
    )
    refseq_fetch_parser.add_argument(
        "--directory",
        default="/refseq/release/complete/",
        help="directory from which files matching '--pattern' argument will be " "downloaded",
    )
    refseq_fetch_parser.add_argument(
        "--pattern", default=r".*protein\.gpff\.gz", help="regex pattern of files to download."
    )
    refseq_fetch_parser.add_argument(
        "--checksum-ftp-path",
        default=r"/refseq/release/release-catalog/.*\.files\.installed",
        help="Path (regex pattern) of a file which contains md5 checksums. Download will "
        "verify that checksum matches for every file if the name is found in the list.",
    )
    refseq_fetch_parser.add_argument(
        "-n",
        "--nr-cpu",
        default=4,
        type=int,
        help="nr of parallel processes to use to download the files from the remote host",
    )
    refseq_fetch_parser.add_argument(
        "-c",
        "--ftp-config",
        help="path to config file for ftp configuration. If not specified, it will try "
        "location in $DARWIN_OMA_RC, or if not set, ~/.omarc. The config files "
        "is not required to exist. The intention is to provide means for a proxy "
        "host for example.",
    )
    refseq_fetch_parser.add_argument(
        "-o",
        "--out",
        default="./",
        help="output directory where the remote files are written to. Defaults to the " "current working directory",
    )

    relevant_taxid_map_parser = subparsers.add_parser(
        "build-taxid-map",
        help="Build taxid mapping to map xrefs with ncbi-taxids information. Includes also subspecies information "
        "and super-species information up to genus level.",
    )
    relevant_taxid_map_parser.set_defaults(func=build_relevant_taxid_mapping)
    relevant_taxid_map_parser.add_argument("--gs-tsv", required=True, help="Path to GS tsv file")
    relevant_taxid_map_parser.add_argument("--db", required=True, help="Path to database in hdf5 format")
    relevant_taxid_map_parser.add_argument("--tax-sqlite", required=False, help="Path to tax-sqlite file")
    relevant_taxid_map_parser.add_argument("--out", required=True, help="Path to output pickle file")

    filter_xref_parser = subparsers.add_parser("filter-xref", help="Filtering xref files")
    filter_xref_parser.set_defaults(func=filter_and_split_xrefs)
    filter_xref_parser.add_argument("--xref", nargs="+", help="Path to input xref files")
    filter_xref_parser.add_argument(
        "--format", required=True, choices=("swiss", "genbank"), help="Format of input xref file"
    )
    filter_xref_parser.add_argument(
        "--out-prefix",
        default="./xref",
        required=False,
        help="Prefix of output xref file. Output files will contain < 30k records, " "all named {prefix}-{03d}.gz",
    )
    filter_xref_parser.add_argument("--tax-map", required=True, help="Path to taxid map file (pickle)")

    map_xref_parser = subparsers.add_parser("map-xref", help="Filtering xref files")
    map_xref_parser.set_defaults(func=map_xrefs)
    map_xref_parser.add_argument("--xref", required=True, help="Path to filtered input xref file")
    map_xref_parser.add_argument(
        "--format", required=True, choices=("swiss", "genbank"), help="Format of input xref file"
    )
    map_xref_parser.add_argument(
        "--source", required=True, choices=("swissprot", "trembl", "refseq"), help="Source of xrefs"
    )
    map_xref_parser.add_argument(
        "--out",
        default="./xref.pkl",
        required=False,
        help="Output file with map results in pickle format",
    )
    map_xref_parser.add_argument("--db", required=True, help="Path to database hdf5 database")
    map_xref_parser.add_argument("--seq-idx-db", required=True, help="Path to sequence index database in hdf5 format")
    map_xref_parser.add_argument("--xref-source-db", required=True, help="Path to xref source hdf5 database")
    map_xref_parser.add_argument("--tax-map", required=True, help="Path to taxid map file (pickle)")

    collect_xref_parser = subparsers.add_parser(
        "collect-xrefs", help="Identify and filter best xrefs matches per source and collect their crossreferences"
    )
    collect_xref_parser.set_defaults(func=collect_xrefs)
    collect_xref_parser.add_argument(
        "--map-results", nargs="+", help="Path to mapped xref pickle files (from map-xref phase)"
    )
    collect_xref_parser.add_argument("--xrefs", nargs="+", help="Path to filtered input xref file")
    collect_xref_parser.add_argument(
        "--format", required=True, choices=("swiss", "genbank"), help="Format of input xref file"
    )
    collect_xref_parser.add_argument(
        "--source", required=True, choices=("swissprot", "trembl", "refseq"), help="Source of xrefs"
    )
    collect_xref_parser.add_argument("--out", required=True, help="Output file with best xref matches per source")

    combine_xref_parser = subparsers.add_parser(
        "combine-xrefs", help="Combine all xref h5 databases into the a single, deduplicated one"
    )
    combine_xref_parser.set_defaults(func=combine_xrefs)
    combine_xref_parser.add_argument("--xrefs", nargs="+", help="Path to input xref files in hdf5 format")
    combine_xref_parser.add_argument("--out", required=True, help="Output path for the combined hdf5 file")

    reduced_xref_parser = subparsers.add_parser(
        "reduced-xrefs", help="Build a reduced set of xrefs for quick search and lookup"
    )
    reduced_xref_parser.set_defaults(func=build_reduced_xrefs)
    reduced_xref_parser.add_argument("--db", required=True, help="Path to database hdf5 database")
    reduced_xref_parser.add_argument("--xrefs", required=True, help="Path to the hdf5 file containing all xrefs")
    reduced_xref_parser.add_argument("--out", required=True, help="Output path for the reduced hdf5 file")
    reduced_xref_parser.add_argument("--nr-procs", type=int, help="Number of processes to use")

    go_import_parser = subparsers.add_parser(
        "import-go", help="Import Gene Ontology ontology (obo file) and annotations (gaf file)"
    )
    go_import_parser.set_defaults(func=import_go)
    go_import_parser.add_argument("--xref-db", required=True, help="Path to xref database in hdf5 format")
    go_import_parser.add_argument("--tax-map", required=True, help="Path to taxid map file (pickle)")
    go_import_parser.add_argument("--obo", required=True, help="Path to input obo file defining gene ontology")
    go_import_parser.add_argument("--gaf", nargs="+", help="Path to one or more gene annotations gaf files")
    go_import_parser.add_argument("--out", required=True, help="Output path for the hdf5 file")

    gen_aux_file_parser = subparsers.add_parser("generate-aux-files", help="Generate auxiliary files")
    gen_aux_file_parser.set_defaults(func=gen_aux_files)
    gen_aux_file_parser.add_argument("--db", required=True, help="Path to database hdf5 database")
    gen_aux_file_parser.add_argument("--out-dir", default="./", help="Path to output directory")

    update_summary_parser = subparsers.add_parser("update-summary", help="Update summary table")
    update_summary_parser.set_defaults(func=store_summary_info)
    update_summary_parser.add_argument(
        "--db", required=True, help="Path to database hdf5 database. This file will be modified!"
    )

    conf = parser.parse_args()
    if not hasattr(conf, "func"):
        parser.print_usage()
        sys.exit(1)
    return conf


def build_database():
    conf = parse_command_line_args()
    logging.basicConfig(
        level=30 - 10 * min(conf.verbose, 2),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    logger.info("Command line options: %s", str(conf))
    if not sys.warnoptions and not getattr(conf, "verbose", 0) >= 1:
        warnings.simplefilter("ignore", category=PerformanceWarning)
        warnings.simplefilter("ignore", category=RuntimeWarning)
    conf.func(conf)


def copy_hdf5_recursive(source_group, target_group, stats=None):
    """
    Recursively merge the contents of source_file into target_file using PyTables.
    Groups with the same name are merged, but datasets (leaves) with the same path are skipped.
    """
    for node in source_group._v_children.values():
        if node._v_name in target_group:
            # Node exists in target, check if it's a group
            if isinstance(node, tables.Group) and isinstance(target_group._v_children[node._v_name], tables.Group):
                logger.info(f"Merging group: {node._v_name}")
                copy_hdf5_recursive(node, target_group._v_children[node._v_name], stats=stats)
            else:
                # Name conflict for non-group items; skipping
                logger.error(f"Skipping duplicate leaf: {node._v_pathname}")
        else:
            # Copy the node (group or leaf)
            logger.info(f"Copying: {node._v_pathname}")
            node._f_copy(target_group, recursive=True, copyuserattrs=True, propindexes=True, stats=stats)


def merge_hdf5():
    parser = ArgumentParser(description="Merge two or more hdf5 files together")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity")
    parser.add_argument("--out", required=True, help="Path to final output hdf5 file")
    parser.add_argument("db", nargs="+", help="Path to an input database in hdf5 format")
    conf = parser.parse_args()
    logging.basicConfig(
        level=30 - 10 * min(conf.verbose, 2),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    logger.info("Command line options: %s", str(conf))
    if not sys.warnoptions and not getattr(conf, "verbose", 0) >= 1:
        warnings.simplefilter("ignore", category=PerformanceWarning)
        warnings.simplefilter("ignore", category=RuntimeWarning)

    if not exists(conf.out):
        out, *rem = sorted(conf.db, key=lambda x: -getsize(x))
        copy(out, conf.out)
        h5_filter = None
    else:
        rem = conf.db
        h5_filter = tables.Filters(complevel=5, complib="blosc2", fletcher32=True)

    stats = {"groups": 0, "leaves": 0, "links": 0, "bytes": 0, "hardlinks": 0}
    with tables.open_file(conf.out, "a", filters=h5_filter) as fout:
        for db in rem:
            with tables.open_file(db, "r") as fin:
                copy_hdf5_recursive(fin.root, fout.root, stats=stats)
    logger.info(stats)
