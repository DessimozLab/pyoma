#!/usr/bin/env python

import csv
import pandas
import os
from tqdm import tqdm
import pyoma.browser.db
import pyoma.browser.hoghelper
import pyoma.hpc
import logging

logger = logging.getLogger("dump-hogids-at-level")


def extract_hogids_at_level(db, level, families, fh_out):
    family_filter = pyoma.browser.hoghelper.HogLevelFilter(db, level)
    for hog_id, level in tqdm(family_filter.analyse_families(families)):
        logger.info("%s %s", hog_id, level)
        fh_out.write(hog_id.decode())
        fh_out.write("\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""Dump all HOG-IDs at a specific taxonomic level. Input comes from an hdf5 file.
                    HOGs that do not reach back as far as the selected reference taxonomic level
                    will be returned as well at their deepest level (if it is a subclade of the
                    selected clade)."""
    )
    parser.add_argument(
        "--out",
        "-o",
        default="hog-dump.txt",
        help="Output file path, defaults to hog-dump.txt",
    )
    parser.add_argument("-v", action="count", default=0, help="Increase verbosity to INFO/DEBUG")
    parser.add_argument("db", help="Path to the hdf5 database file")
    parser.add_argument("level", help="Level at which to produce the groups")
    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10 * min(conf.v, 2))

    db = pyoma.browser.db.Database(conf.db)
    families = range(db.get_nr_toplevel_hogs() + 1)
    with open(conf.out, "wt", encoding="utf-8") as fh_out:
        extract_hogids_at_level(db, conf.level, families, fh_out)
