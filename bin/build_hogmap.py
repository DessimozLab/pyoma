import argparse
import logging
import pyoma.browser.hogidmap

logger = logging.getLogger("build-hog-map")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="build database to find related hogs in new version")
    parser.add_argument("-v", action="count", default=0)
    subparser = parser.add_subparsers(required=True)

    parser_lsh = subparser.add_parser("LSH")
    parser_lsh.add_argument(
        "--db",
        required=True,
        help="Path to database for which the LSH of HOGs should be build",
    )
    parser_lsh.add_argument("--out", required=True, help="Path where the lsh should be written to")
    parser_lsh.add_argument(
        "--nr-procs",
        type=int,
        help="Nr of processes to use. defaults to all available cores of the node",
    )
    parser_lsh.set_defaults(stage="build_lsh")
    parser_map = subparser.add_parser("map")
    parser_map.add_argument("--target", required=True, help="Path to the LSH database of the target release")
    parser_map.add_argument("--old", nargs="+", help="Path to the LSH databases of the old releases.")
    parser_map.add_argument(
        "--out",
        required=True,
        help="Path where the hogmap database should be written to",
    )
    parser_map.set_defaults(stage="map")

    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10 * min(conf.v, 2))
    logger.info(conf)

    if conf.stage == "build_lsh":
        pyoma.browser.hogidmap.compute_minhashes_for_db(conf.db, conf.out, nr_procs=conf.nr_procs)
    elif conf.stage == "map":
        pyoma.browser.hogidmap.compare_versions(conf.out, conf.target, *conf.old)

    print("done... bye bye")
