#!/usr/bin/env python

import pyoma.browser.db
import Bio.SeqIO
import logging
import csv
logger = logging.getLogger('hog_mapper')



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Map sequences to HOGs")
    parser.add_argument('hdf5', help="Path to the hdf5 database")
    parser.add_argument('fasta', nargs="+", help="File(s) with fasta formatted protein sequences")
    parser.add_argument('--out', default="map2hog", help="filename prefix of output file.")
    parser.add_argument('-n', '--nr_proc', help="Nr of processes to use")
    parser.add_argument('-v', action='count', default=0, help="Increase verbosity to INFO/DEBUG")
    conf = parser.parse_args()
    logging.basicConfig(level=30 - min(conf.v, 2))

    db = pyoma.browser.db.Database(conf.hdf5)
    hog_mapper = pyoma.browser.db.SimpleSeqToHOGMapper(db)
    outfn = conf.out + ".txt"
    with open(outfn, 'wt') as fout:
        csv_writer = csv.writer(fout, delimiter="\t")
        for infn in conf.fasta:
            with open(infn, 'rt') as fin:
                seqs = Bio.SeqIO.parse(fin, 'fasta')
                it = hog_mapper.imap_sequences(seqs)

                for map_res in it:
                    if map_res.target.hog_family_nr != 0:
                        hog_id = db.format_hogid(map_res.target.hog_family_nr)
                    else:
                        hog_id = "n/a"
                    csv_writer.writerow([map_res.query, map_res.target.omaid,
                                         map_res.target.entry_nr == map_res.closest_entry_nr,
                                         hog_id, map_res.distance, map_res.score])
