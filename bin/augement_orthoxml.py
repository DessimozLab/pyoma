#!/usr/bin/env python

from pyoma.browser.convert import HogConverter
import tables
import logging

logger = logging.getLogger("augment_orthoxml")


def augment_orthoxml(h5, taxonomy, orthoxml, release, out=None):
    if out is None:
        out = orthoxml + ".augmented"
    entry_tab = h5.get_node("/Protein/Entries")
    tax_tab = h5.get_node("/Taxonomy")
    tax_2_code = {int(row["NCBITaxonId"]): row["UniProtSpeciesCode"].decode() for row in h5.get_node("/Genome")}
    hog_converter = HogConverter(entry_tab, release, tax_tab, tax_2_code)
    hog_converter.attach_newick_taxonomy(taxonomy)
    levels = hog_converter.convert_file(orthoxml, store=out)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Augment existing orthoxml file")
    parser.add_argument("--hdf5", required=True, help="Path to the hdf5 database")
    parser.add_argument("--orthoxml-file", required=True, help="input orthoxml file")
    parser.add_argument("--out", help="path to output file. defaults to inputfile + '.augmented'")
    parser.add_argument("--taxonomy", required=True, help="Taxonomy file in newick format")
    parser.add_argument(
        "--release",
        "-r",
        required=True,
        help="Release number as a single character (for hogs)",
    )
    conf = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    logger.info("Params: %s", conf)

    with tables.open_file(conf.hdf5, "r") as h5:
        augment_orthoxml(h5, conf.taxonomy, conf.orthoxml_file, conf.release, conf.out)
