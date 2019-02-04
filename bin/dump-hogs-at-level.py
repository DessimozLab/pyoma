#!/usr/bin/env python

import csv
import pandas
import os
from tqdm import tqdm
import pyoma.browser.db
import pyoma.browser.hoghelper
import pyoma.hpc
import logging
logger = logging.getLogger('dump-hogs-at-level')


def load_xref(db, idtype=None):
    if not isinstance(db, pyoma.browser.db.Database):
        raise TypeError("db must be a hdf5 oma database")

    xref_tab = db.db.get_node('/XRef')
    xrefs = pandas.DataFrame(xref_tab.read_where(
        '(XRefSource == {})'.format(xref_tab.get_enum('XRefSource')[idtype])))
    xrefs.drop(columns=['XRefSource', 'Verification'], inplace=True)
    xrefs['XRefId'] = xrefs['XRefId'].str.decode("utf-8")
    xrefs.rename(columns={'XRefId': idtype}, inplace=True)
    return xrefs


def extract_hogs_at_level(db, level, families, fh_out, xref_types=None):
    family_filter = pyoma.browser.hoghelper.HogLevelFilter(db, level)

    cols = ["HOG", "Level", "OmaID"]
    if xref_types is not None:
        cols.extend(xref_types)
        for k, xtype in enumerate(xref_types):
            xref_df = load_xref(db, xtype)
            logger.info('loaded {} xrefs for {}'.format(len(xref_df), xtype))
            if k == 0:
                xrefs_df = xref_df
            else:
                no_dups = xref_df.drop_duplicates(keep='first', subset='EntryNr')
                logger.info('dropped {} of {} rows that had several {} xrefs per EntryNr'
                            .format(len(xref_df)-len(no_dups), len(xref_df), xtype))
                xrefs_df = xrefs_df.merge(no_dups, how='outer', on='EntryNr')


    csv_writer = csv.DictWriter(fh_out, cols, extrasaction="ignore", delimiter="\t")
    csv_writer.writeheader()
    for hog_id, level in tqdm(family_filter.analyse_families(families)):
        logger.info("{} {}".format(hog_id, level))
        df = pandas.DataFrame(db.member_of_hog_id(hog_id, level))
        df['OmaID'] = df['EntryNr'].apply(db.id_mapper['Oma'].map_entry_nr)
        df['HOG'] = hog_id.decode() if isinstance(hog_id, bytes) else hog_id
        df['Level'] = level.decode() if isinstance(level, bytes) else level
        if xref_types is not None:
            df = df.merge(xrefs_df, how='left', on="EntryNr")
            #df.rename(columns={'XRefId': xref_type}, inplace=True)
        csv_writer.writerows(df.to_dict('records'))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Dump all HOGs at a specific taxonomic level. Input comes from an hdf5 file.
                    HOGs that do not reach back as far as the selected reference taxonomic level 
                    will be returned as well at their deepest level (if it is a subclade of the 
                    selected clade).""")
    parser.add_argument('--xref-type', nargs='*',
                        help="CrossReference ID used in output. The value must be matching one of the "
                             "enum values in the XRef table, e.g. SourceAC")
    parser.add_argument('--out', '-o', default="hog-dump.txt",
                        help="Output file path, defaults to hog-dump.txt")
    parser.add_argument('-n', '--nr_proc', type=int, help="Nr of processes to use")
    parser.add_argument('-p', '--procnr', type=int,
                        help="This process nr. If not passed explictly, it tries to figure out from jobarray info")
    parser.add_argument('-v', action='count', default=0, help="Increase verbosity to INFO/DEBUG")
    parser.add_argument('db', help="Path to the hdf5 database file")
    parser.add_argument('level', help="Level at which to produce the groups")
    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10 * min(conf.v, 2))
    nr_procs = conf.nr_proc
    if nr_procs is None:
        nr_procs = int(os.getenv('NR_PROCESSES', "1"))
    pInf = pyoma.hpc.detect_hpc_jobarray(nr_procs, this_proc_nr=conf.procnr)
    outfn = pInf.modify_filename(conf.out)
    logger.info(pInf)

    db = pyoma.browser.db.Database(conf.db)
    families = filter(lambda x: pInf.is_my_job(x), range(db.get_nr_toplevel_hogs()))
    with open(outfn, 'w') as fh_out:
        extract_hogs_at_level(db, conf.level, families, fh_out, xref_types=conf.xref_type)
