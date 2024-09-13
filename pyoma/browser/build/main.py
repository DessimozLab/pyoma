import logging
import sys
import warnings

import tables
from tables import PerformanceWarning
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# from .. import convert
from .builder import DBBuilder, OmaGroupsProvider, XrefStorer
from . import hogconvert

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
        hogconvert.PerFamilyHOGObserver(parser, hog_h5.store_orthoxml, hog_h5.store_orthoxml_augmented)
        hogconvert.parse_orthoxml(conf.orthoxml, parser)


def phase_vps(conf):
    with tables.open_file(conf.db, "r") as db:
        genomes = db.get_node("/Genome").read()
    with DBBuilder(conf.hdf5_out, mode="write", logger=logger, complib="blosc") as out:
        out.add_orthologs(conf.vps_base, genomes=genomes)


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
