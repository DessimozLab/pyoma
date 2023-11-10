import argparse
import tables
import pyoma.browser.xref_contrib


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="build reduced xref lookup data structure")
    parser.add_argument("--db", required=True, help="Path to the hdf5 database")
    parser.add_argument(
        "--out",
        default="reduced_xrefs.h5",
        help="Path output hdf5, must be different from --db",
    )
    parser.add_argument(
        "--nr-procs",
        type=int,
        help="Nr of processes to use. defaults to all available cores of the node",
    )
    conf = parser.parse_args()
    print(conf)
    if conf.db == conf.out:
        raise Exception("db and out parameters must be different")

    pyoma.browser.xref_contrib.reduce_xrefs(conf.db, conf.out, nr_procs=conf.nr_procs)
    print("done... bye bye")
