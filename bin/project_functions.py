#!/usr/bin/env python
import os

import pyoma.browser.db
from pyoma.common import auto_open
import pyoma.hpc
import Bio.SeqIO
import logging
import sys
import signal

logger = logging.getLogger("mapper")


def handler(signum, frame):
    logger.info(
        "signal handler called with signal {}. Raise a KeybordInterupt Exception".format(
            signum
        )
    )
    raise KeyboardInterrupt("time is up")


signal.signal(signal.SIGUSR2, handler)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Map sequences to HOGs")
    parser.add_argument("hdf5", help="Path to the hdf5 database")
    parser.add_argument(
        "fasta",
        nargs="+",
        help="File(s) with fasta formatted protein sequences. Checkpointing does only work with a single input file so far.",
    )
    parser.add_argument(
        "--out",
        default="function.gaf",
        help="filename output file. Will be modified if run with more than one process",
    )
    parser.add_argument("-n", "--nr_proc", type=int, help="Nr of processes to use")
    parser.add_argument("-p", "--procnr", type=int, help="This process nr")
    parser.add_argument(
        "-v", action="count", default=0, help="Increase verbosity to INFO/DEBUG"
    )
    parser.add_argument(
        "-t",
        "--target",
        default=None,
        help="Target species. If not set, no restriction applies.",
    )
    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10 * min(conf.v, 2))
    nr_procs = conf.nr_proc
    if nr_procs is None:
        nr_procs = int(os.getenv("NR_PROCESSES", "1"))
    pInf = pyoma.hpc.detect_hpc_jobarray(nr_procs, this_proc_nr=conf.procnr)
    logger.info(pInf)
    outfn = pInf.modify_filename(conf.out)

    db = pyoma.browser.db.Database(conf.hdf5)
    mapper = pyoma.browser.db.FastMapper(db)

    last_id = None
    if os.path.exists(outfn + ".ckpt"):
        with open(outfn + ".ckpt") as fh:
            last_id = fh.read()

    with auto_open(outfn, "at") as fout:
        for infn in conf.fasta:
            with auto_open(infn, "rt") as fin:
                seqs = Bio.SeqIO.parse(fin, "fasta")
                myjobs = filter(lambda srec: pInf.is_my_job(srec.id), seqs)
                # forward in case we have a
                if last_id is not None:
                    for xx in myjobs:
                        if xx.id == last_id:
                            break
                try:
                    mapper.write_annotations(fout, myjobs, target_species=conf.target)

                except KeyboardInterrupt as e:
                    logger.info("received KeyboardInterrupt: {}".format(e))
                    logger.info("writing checkpoint")
                    with open(outfn + ".ckpt", "w") as ckpt:
                        ckpt.write(last_id)
                    logger.info("quiting with 99")
                    if len(conf.fasta) > 1:
                        logger.warning(
                            "checkpointing does only properly work with a single input file. You used more than one file. please investigate how to recover"
                        )
                        sys.exit(98)
                    sys.exit(99)
