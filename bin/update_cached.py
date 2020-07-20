import pyoma.browser.compute_cache
import logging
import argparse

logger = logging.getLogger("cache-update")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Update cached/precomputed data for browser"
    )
    parser.add_argument("hdf5", help="Path to the hdf5 database")
    parser.add_argument("-v", default=0, action="count", help="increase verbosity")
    parser.add_argument(
        "-n", "--nr-procs", default=None, type=int, help="nr of processes to use"
    )
    parser.add_argument(
        "-c", "--cache-file", default=None, help="Path to a temporary cache file"
    )
    parser.add_argument(
        "-f",
        "--overwrite",
        action="store_true",
        default=False,
        help="Overwrite existing cache values",
    )
    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10 * min(conf.v, 2))

    cache = pyoma.browser.compute_cache.compute_and_store_cached_data(
        conf.hdf5,
        "/Protein/OrthologsCountCache",
        conf.nr_procs,
        force=conf.overwrite,
        tmp_cache=conf.cache_file,
    )
