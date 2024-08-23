import logging
import sys
import warnings
from tables import PerformanceWarning
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .. import convert
from .builder import DBBuilder, OmaGroupsProvider, XrefStorer

logger = logging.getLogger(__name__)


def phase_genomes(conf):
    with DBBuilder(conf.db, mode="write", logger=logger) as db:
        db.add_version(conf.rel_char)
        db.add_taxonomy(conf.tax_tsv)
        db.add_species_data(conf.gs_tsv)
        with XrefStorer(conf.xref_db) as xref_storer:
            db.add_proteins(conf.genomes, OmaGroupsProvider(conf.oma_groups), xref_collector=xref_storer)


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
