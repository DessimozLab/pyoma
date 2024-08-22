import logging
import sys
import warnings
from tables import PerformanceWarning
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .. import convert
from .builder import DBBuilder

logger = logging.getLogger(__name__)


def phase_genomes(conf):
    with DBBuilder(conf.db, mode="write", logger=logger) as db:
        db.add_species_data(conf.gs_tsv, conf.tax_tsv)


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
    genomes_parser.add_argument(
        "--genomes", required=True, nargs="+", help="List of genome files (json) containing essential data"
    )

    conf = parser.parse_args()
    logging.basicConfig(
        level=30 - 10 * min(conf.verbose, 2),
        format="%(asctime)s %(levelname)s %(name): %(message)s",
    )
    logger.info("Command line options: %s", str(conf))
    if hasattr(conf, "func"):
        if not sys.warnoptions and not getattr(conf, "verbose", 0) >= 1:
            warnings.simplefilter("ignore", category=PerformanceWarning)
            warnings.simplefilter("ignore", category=RuntimeWarning)
        conf.func(conf)
    else:
        parser.print_usage()


def build_database():
    conf = parse_command_line_args()
