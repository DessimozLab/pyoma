#!/usr/bin/env python3

from pyoma.browser.hogprofile import compute_profiles
import logging
import argparse

logger = logging.getLogger("compute-profiles")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute hog profiles for browser")
    parser.add_argument("hdf5", help="Path to the hdf5 database")
    parser.add_argument("-v", default=0, action="count", help="increase verbosity")
    parser.add_argument(
        "-n", "--nr-procs", default=None, type=int, help="nr of processes to use"
    )
    parser.add_argument(
        "-m",
        "--min-hogsize",
        type=int,
        default=70,
        help="minimal nr of species in a hog to be considered",
    )
    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10 * min(conf.v, 2))

    cache = compute_profiles(
        conf.hdf5, min_hogsize=conf.min_hogsize, nr_procs=conf.nr_procs
    )
