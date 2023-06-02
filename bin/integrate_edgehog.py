#!/usr/bin/env python3
import tables
import logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Integrate results from EdgeHOG into oma browser database"
    )
    parser.add_argument("-v", action="count", default=0, help="increase verbosity")
    parser.add_argument("--omadb", required=True, help="path to OmaServer.h5 file")
    parser.add_argument(
        "--edgehog", required=True, help="Path to edgehog result in hdf5 format"
    )

    conf = parser.parse_args()
    log_level = 30 - (10 * min(conf.v, 2))
    logging.basicConfig(level=log_level)
    logger.info("Params: {}".format(conf))
    with tables.open_file(conf.edgehog, "r") as edge_h5, tables.open_file(
        conf.omadb, "a"
    ) as h5:
        for node in edge_h5.walk_nodes("/", tables.Table):
            path, name = node._v_pathname.rsplit("/", 1)
            target_node = h5.get_node(path)
            logger.debug(f"copy {name} from {path} to {target_node}")
            node._f_copy(target_node, name)
