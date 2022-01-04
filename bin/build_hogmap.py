import pyoma.browser.hogidmap
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="find related hogs in new version")
    parser.add_argument("--target", required=True, help="Path to the target database")
    parser.add_argument(
        "--old", nargs="+", help="Path to databases that should be mapped"
    )
    parser.add_argument("--nr-procs", type=int, 
                        help="Nr of processes to use. defaults to all available cores of the node")
    conf = parser.parse_args()
    print(conf)

    pyoma.browser.hogidmap.build_lookup(conf.target, conf.old, nr_procs=conf.nr_procs)
    print("done... bye bye")
