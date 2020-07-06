#!/usr/bin/env python

from pyoma.browser.convert import HogConverter
import tables
import logging

logger = logging.getLogger("augment_orthoxml")


def augment_orthoxml(h5, taxonomy, orthoxml):
    entry_tab = h5.get_node("/Protein/Entries")
    tax_tab = h5.get_node("/Taxonomy")
    hog_converter = HogConverter(entry_tab, tax_tab)
    hog_converter.attach_newick_taxonomy(taxonomy)
    out_orthoxml = orthoxml + ".augmented"
    levels = hog_converter.convert_file(orthoxml, store=out_orthoxml)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Augment existing orthoxml file")
    parser.add_argument("--hdf5", required=True, help="Path to the hdf5 database")
    parser.add_argument("--orthoxml-file", required=True, help="input orthoxml file")
    parser.add_argument(
        "--taxonomy", required=True, help="Taxonomy file in newick format"
    )
    conf = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    logger.info("Params: {}".format(conf))

    with tables.open_file(conf.hdf5, "r") as h5:
        augment_orthoxml(h5, conf.taxonomy, conf.orthoxml_file)
