#!/usr/bin/env python
from distutils.spawn import find_executable
import pyoma.browser.convert
import pyoma.browser.convert_omastandalone
import argparse
import os
import sys


def main(args):
    # Check that darwin is installed on the system.
    if not find_executable('darwin'):
        raise RuntimeError('oma2hdf requires darwin: can\'t find darwin on your system.')

    # Argument parsing
    parser = argparse.ArgumentParser("Convert a OMA Browser release into hdf5 format.")
    parser.add_argument('-r', '--release',
                        help="path to release. If not set, DARWIN_BROWSERDATA_PATH "
                             "from the environment is used.")
    parser.add_argument('out', default='OmaServer.h5',
                        help="name of the hdf5 database that is created. The file is "
                             "stored in the same path as the release is stored, i.e. "
                             "see option -r/--release and it's default.")
    parser.add_argument('-s', '--standalone', action='store_true',
                        help="a flag which needs to be set if you intend to import "
                             "the results from an OmaStandalone run.")
    parser.add_argument('--no-domains', action='store_false',
                        help="do not include CATH domain information. This flag is only "
                             "considered for oma standalone imports. (see -s/--standalone flag")

    options = parser.parse_args(args)
    if options.standalone:
        pyoma.browser.convert_omastandalone.import_oma_run(
            options.release, options.out, add_domains=options.no_domains)
    else:
        if options.release:
            os.environ['DARWIN_BROWSERDATA_PATH'] = options.release
        print(options.out)

        pyoma.browser.convert.main(options.out)


if __name__ == '__main__':
    main(sys.argv[1:])
