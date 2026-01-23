import logging
import sys
import warnings
from contextlib import ExitStack
from collections import defaultdict
from textwrap import dedent

import tables
from tables import PerformanceWarning
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from .main import setup_logging
from ..db import Database
from ..build.taxonomy_alignment import produce_mapping_and_download_files
from ...common import auto_open

logger = logging.getLogger(__name__)


def open_xref_outputs(args):
    """Return dict: SourceXRef -> open file handle"""
    # Map sources to paths
    source_to_path = {
        "UniProtKB/SwissProt": args.out_uniprot,
        "UniProtKB/TrEMBL": args.out_uniprot,
        "Ensembl Protein": args.out_ensembl,
        "Ensembl Gene": args.out_ensembl,
        "Ensembl Transcript": args.out_ensembl,
        "RefSeq": args.out_refseq,
        "NCBI": args.out_ncbi,
        "EntrezGene": args.out_entrez,
    }

    # Open each unique path only once
    unique_paths = {p for p in source_to_path.values() if p}
    stack = ExitStack()
    path_to_handle = {p: stack.enter_context(auto_open(p, "wt")) for p in unique_paths}

    # Map each source to the corresponding handle
    outputs = {source: path_to_handle[path] for source, path in source_to_path.items() if path}
    return outputs, stack


def dump_xrefs(conf):
    with Database(conf.db) as db:
        if conf.xref_db is not None:
            xref_db = tables.open_file(conf.xref_db)
            xref_tab = xref_db.get_node("/XRef")
        else:
            xref_tab = db.db.get_node("/XRef")
        oma_id_mapper = db.id_mapper["OMA"]
        source_enum = xref_tab.get_enum("XRefSource")

        outputs, stack = open_xref_outputs(conf)
        with stack:
            for row in xref_tab.iterrows():
                source = source_enum(row["XRefSource"])
                if source not in outputs:
                    continue
                omaid = oma_id_mapper.map_entry_nr(row["EntryNr"])
                outputs[source].write(f"{omaid}\t{row['XRefId'].decode('utf-8')}\n")

        if conf.xref_db is not None:
            xref_db.close()


def dump_go(conf):
    with Database(conf.db) as db:
        oma_id_mapper = db.id_mapper["OMA"]
        with auto_open(conf.out_go, "wt") as go_fh:
            STEP = 100_000
            go_fh.write(
                dedent(
                    f"""\
            # Mapping of OMA IDs to Gene Ontology (GO) annotations for OMA release {db.get_release_name()}
            # Format: OMA ID<tab>GO term<tab>Evidence code<tab>Reference codes (comma-separated)
            """
                )
            )
            for enr in range(1, db.nr_proteins + 1, STEP):
                agg = defaultdict(set)
                for row in db.get_gene_ontology_annotations(enr, stop=enr + STEP):
                    key = (
                        oma_id_mapper.map_entry_nr(row["EntryNr"]),
                        f"GO:{row['TermNr']:07d}",
                        row["Evidence"].decode("utf-8"),
                    )
                    agg[key].add(row["Reference"].decode("utf-8"))
                for (omaid, term, evi), refs in agg.items():
                    refs_str = ",".join(refs)
                    go_fh.write(f"{omaid}\t{term}\t{evi}\t{refs_str}\n")


def dump_species_and_taxmapping(conf):
    produce_mapping_and_download_files(conf.db, conf.tax_sqlite, conf.out_species, conf.out_tax_mapping)


def parse_command_line_args():
    parser = ArgumentParser(description="Dump various files from the OMA database")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity")

    subparsers = parser.add_subparsers(title="Commands")

    xrefs_parser = subparsers.add_parser(
        "xrefs", help="Dump cross-reference mappings", formatter_class=ArgumentDefaultsHelpFormatter
    )
    xrefs_parser.set_defaults(func=dump_xrefs)
    xrefs_parser.add_argument("--db", required=True, help="Path to database")
    xrefs_parser.add_argument(
        "--xref-db", required=False, help="Path to xref database. If not provided, the main database is used."
    )
    xrefs_parser.add_argument("--out-uniprot", help="Output path for the uniprot cross-reference mappings (if enabled)")
    xrefs_parser.add_argument("--out-ensembl", help="Output path for the ensembl cross-reference mappings (if enabled)")
    xrefs_parser.add_argument("--out-refseq", help="Output path for the refseq cross-reference mappings (if enabled)")
    xrefs_parser.add_argument("--out-ncbi", help="Output path for the ncbi cross-reference mappings (if enabled)")
    xrefs_parser.add_argument("--out-entrez", help="Output path for the entrez cross-reference mappings (if enabled)")

    go_parser = subparsers.add_parser("go", help="Dump GO annotations", formatter_class=ArgumentDefaultsHelpFormatter)
    go_parser.set_defaults(func=dump_go)
    go_parser.add_argument("--db", required=True, help="Path to database")
    go_parser.add_argument("--out-go", help="Output path for the GO annotations")

    species_parser = subparsers.add_parser(
        "species", help="Dump species information", formatter_class=ArgumentDefaultsHelpFormatter
    )
    species_parser.set_defaults(func=dump_species_and_taxmapping)
    species_parser.add_argument("--db", required=True, help="Path to database")
    species_parser.add_argument("--tax-sqlite", required=True, help="Path to the taxonomy.sqlite file")
    species_parser.add_argument(
        "--out-species", required=True, help="Output path for the species information in tsv format"
    )
    species_parser.add_argument(
        "--out-tax-mapping", required=True, help="Path where the GTDB / NCBI taxonomy mapping is stored as pickle file."
    )

    conf = parser.parse_args()
    if not hasattr(conf, "func"):
        parser.print_usage()
        sys.exit(1)
    return conf


def dump_files():
    conf = parse_command_line_args()
    setup_logging(conf)
    logger.info("Command line options: %s", str(conf))
    if not sys.warnoptions and not getattr(conf, "verbose", 0) >= 1:
        warnings.simplefilter("ignore", category=PerformanceWarning)
        warnings.simplefilter("ignore", category=RuntimeWarning)
    conf.func(conf)


if __name__ == "__main__":
    dump_files()
