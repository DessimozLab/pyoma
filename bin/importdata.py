#!/usr/bin/env python
from distutils.spawn import find_executable
import pyoma.browser.convert
import pyoma.browser.convert_omastandalone
import argparse
import os
import sys


def main(args):
    # Check that darwin is installed on the system.
    if not find_executable("darwin"):
        raise RuntimeError("oma2hdf requires darwin: can't find darwin on your system.")

    # Argument parsing
    parser = argparse.ArgumentParser("Convert a OMA Browser release into hdf5 format.")
    parser.add_argument(
        "-r",
        "--release",
        help="path to release. If not set, DARWIN_BROWSERDATA_PATH "
        "from the environment is used.",
    )
    parser.add_argument(
        "out",
        default="OmaServer.h5",
        help="name of the hdf5 database that is created. The file is "
        "stored in the same path as the release is stored, i.e. "
        "see option -r/--release and it's default.",
    )
    parser.add_argument(
        "-s",
        "--standalone",
        action="store_true",
        help="a flag which needs to be set if you intend to import "
        "the results from an OmaStandalone run.",
    )
    parser.add_argument(
        "-d",
        "--domains",
        nargs="+",
        help="absolute path or url to domain annotations in mdas.csv format. "
        "Not specifying any domains is equivalent to the --no-domains "
        "option.",
    )
    parser.add_argument(
        "--hog-release-char",
        required=("--standalone" not in sys.argv and "-s" not in sys.argv),
        help="A single character indicating the release in the hog ids. "
        "This argument is ignored for oma standalone imports",
    )
    parser.add_argument(
        "-v",
        default=0,
        action="count",
        help="Increase verbosity level to INFO or DEBUG level",
    )

    options = parser.parse_args(args)
    log_level = 30 - (10 * min(options.v, 2))
    if options.standalone:
        pyoma.browser.convert_omastandalone.import_oma_run(
            options.release, options.out, domains=options.domains, log_level=log_level
        )
    else:
        if options.release:
            os.environ["DARWIN_BROWSERDATA_PATH"] = options.release
        print(options.out)
        pyoma.browser.convert.main(
            options.out, domains=options.domains, log_level=log_level
        )


if __name__ == "__main__":
    main(sys.argv[1:])
