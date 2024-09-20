import collections
import itertools
import logging
import sys
import warnings
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import ete3
import pandas
import tables
from tables import PerformanceWarning
import omataxonomy


from .builder import DBBuilder, OmaGroupsProvider, XrefStorer
from . import hogconvert
from pyoma.browser.build import xref as xref_build
from ..convert import (
    iter_domains,
    filter_duplicated_domains,
    only_pfam_or_cath_domains,
    CathDomainNameParser,
    PfamDomainNameParser,
)

logger = logging.getLogger(__name__)


def phase_genomes(conf):
    with DBBuilder(conf.db, mode="write", logger=logger) as db:
        db.add_version(conf.rel_char)
        db.add_taxonomy(conf.tax_tsv)
        db.add_species_data(conf.gs_tsv)
        with XrefStorer(conf.xref_db) as xref_storer:
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


def fetch_refseq(conf):
    xref_build.fetch(**vars(conf))


def filter_and_split_xrefs(conf):
    gs = pandas.read_csv(conf.gs_tsv, sep="\t")
    ncbi_taxids = set(gs["OriginalNCBITaxonId"])
    relevant_taxids = xref_build.load_relevant_taxids(ncbi_taxids, omataxonomy.Taxonomy(conf.tax_sqlite))
    with xref_build.ChunkWriter(conf.out_prefix) as writer:
        xref_build.filter_records_on_taxids(conf.xref, writer, relevant_taxids, conf.format)


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
    vp_parser.add_argument("--vps-base", required=True, help="Folder where all the pairwise orthologs are stored")
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

    filter_xref_parser = subparsers.add_parser("filter-xref", help="Filtering xref files")
    filter_xref_parser.set_defaults(func=filter_and_split_xrefs)
    filter_xref_parser.add_argument("--xref", required=True, help="Path to input xref file")
    filter_xref_parser.add_argument(
        "--format", required=True, choices=("swiss", "genbank"), help="Format of input xref file"
    )
    filter_xref_parser.add_argument(
        "--out-prefix",
        default="./xref",
        required=False,
        help="Prefix of output xref file. Output files will contain < 30k records, " "all named {prefix}-{03d}.gz",
    )
    filter_xref_parser.add_argument("--gs-tsv", required=True, help="Path to GS tsv file")
    filter_xref_parser.add_argument("--tax-sqlite", required=False, help="Path to tax-sqlite file")

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
