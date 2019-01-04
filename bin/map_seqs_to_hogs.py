#!/usr/bin/env python
import os

import pyoma.browser.db
import pyoma.hpc
import Bio.SeqIO
import logging
import csv
import hashlib
logger = logging.getLogger('hog_mapper')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Map sequences to HOGs")
    parser.add_argument('hdf5', help="Path to the hdf5 database")
    parser.add_argument('fasta', nargs="+", help="File(s) with fasta formatted protein sequences")
    parser.add_argument('--out', default="map2hog.txt",
                        help="filename output file. Will be modified if run with more than one process")
    parser.add_argument('-n', '--nr_proc', type=int, help="Nr of processes to use")
    parser.add_argument('-p', '--procnr', type=int, help="This process nr")
    parser.add_argument('-v', action='count', default=0, help="Increase verbosity to INFO/DEBUG")
    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10*min(conf.v, 2))
    nr_procs = conf.nr_proc
    if nr_procs is None:
        nr_procs = int(os.getenv('NR_PROCESSES', "1"))
    pInf = pyoma.hpc.detect_hpc_jobarray(nr_procs, this_proc_nr=conf.procnr)
    logger.info(pInf)
    outfn = pInf.modify_filename(conf.out)

    db = pyoma.browser.db.Database(conf.hdf5)
    hog_mapper = pyoma.browser.db.SimpleSeqToHOGMapper(db)

    with open(outfn, 'wt') as fout:
        csv_writer = csv.writer(fout, delimiter="\t")
        if pInf.this_proc_nr == 1:
            csv_writer.writerow(['query', 'target', 'is_main_isoform', 'HOG', 'PAM_distance', 'Alignment_score'])
            
        for infn in conf.fasta:
            with open(infn, 'rt') as fin:
                seqs = Bio.SeqIO.parse(fin, 'fasta')
                myjobs = filter(lambda srec: pInf.is_my_job(srec.id), seqs)
                it = hog_mapper.imap_sequences(myjobs)

                for map_res in it:
                    if map_res.target.hog_family_nr != 0:
                        hog_id = db.format_hogid(map_res.target.hog_family_nr)
                    else:
                        hog_id = "n/a"
                    csv_writer.writerow([map_res.query, map_res.target.omaid,
                                         map_res.target.entry_nr == map_res.closest_entry_nr,
                                         hog_id, map_res.distance, map_res.score])
