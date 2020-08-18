import pyoma.browser.ancestral_synteny
import logging
import argparse

logger = logging.getLogger("ancestral_synteny")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute ancestral synteny")
    parser.add_argument("--hdf5", required=True, help="Path to the hdf5 database")
    parser.add_argument("--orthoxml", required=True, help="Path to the orthoxml file")
    parser.add_argument(
        "--tree", required=True, help="Path to the species tree in newick format"
    )
    parser.add_argument("-v", default=0, action="count", help="increase verbosity")
    parser.add_argument(
        "-n", "--nr-procs", default=None, type=int, help="nr of processes to use"
    )
    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10 * min(conf.v, 2))

    res = pyoma.browser.ancestral_synteny.infer_synteny(
        conf.orthoxml, conf.hdf5, conf.tree
    )
