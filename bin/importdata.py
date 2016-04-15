#!/usr/bin/env python

import pyoma.browser.convert
import pyoma.browser.convert_omastandalone
import argparse
import os
import sys

def main(args):
    parser = argparse.ArgumentParser("Convert a OMA Browser release into hdf5 format.")
    parser.add_argument('-r', '--release', help="path to release. If not set, DARWIN_BROWSERDATA_PATH"
                "from the environment is used.")
    parser.add_argument('out', default='OmaServer.h5')
    parser.add_argument('-v', '--verbose', action='count')
    parser.add_argument('-s', '--standalone', action='store_true',
                        help="to be set if import from an OmaStandalone run")

    options = parser.parse_args(args)
    if options.standalone:
        pyoma.browser.convert_omastandalone.import_oma_run(options.release, options.out)
    else:
        if options.release:
            os.environ['DARWIN_BROWSERDATA_PATH'] = options.release
        print(options.out)

        pyoma.browser.convert.main(options.out)

if __name__ == '__main__':
    main(sys.argv[1:])
