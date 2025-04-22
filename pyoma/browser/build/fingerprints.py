from typing import Optional, List, Dict
import logging
from collections import Counter
import numpy
from tqdm import tqdm

from ..db import Database, SequenceSearch, read_table_where

logger = logging.getLogger(__name__)


def kmers(seqs: List, k: int) -> Counter:
    """
    Count the kmers in a set of sequences.
    :param seqs:
    :param k:
    :return:
    """
    kmer_counts = Counter()
    for seq in seqs:
        for i in range(len(seq) - k - 1):  # Ensure at least k chars remain
            kmer = seq[i : i + k]
            kmer_counts[kmer] += 1
    return kmer_counts


def find_fingerprints(
    db_path: str, suffix_path: Optional[str] = None, ogs: Optional[List[int]] = None
) -> Dict[int, str]:
    """
    Find fingerprints in a given set of oma groups.

    :param db_path: path to the OMA database
    :param suffix_path: path to the suffix array file
    :param List[int] ogs: list of oma groups to process. if None, all groups are processed
    :return: dictionary with the oma groups as keys and the fingerprints as values
    """
    db = Database(db_path)
    searcher = SequenceSearch(db, seq_idx_fpath=suffix_path)
    # read the sequence buffer into memory
    searcher.seq_buff = searcher.seq_buff[:]
    pe_tab = db.db.get_node("/Protein/Entries")

    og_iter = range(1, db.get_nr_oma_groups() + 1) if ogs is None else ogs
    fingerprints = {}
    for og in tqdm(og_iter, desc="Finding fingerprints"):
        og_entries = read_table_where(pe_tab, "(OmaGroup == og)", condvars={"og": og})
        seqs = [db.get_sequence(e) for e in og_entries]
        kmers_cnt = kmers(seqs, 7)

        for km, cnt in kmers_cnt.most_common():
            enrs = searcher.exact_search(km, only_full_length=False, is_sanitised=True)
            if numpy.isin(enrs, og_entries["EntryNr"]).all():
                fingerprints[og] = km.decode()
                break
        else:
            fingerprints[og] = "n/a"
    return fingerprints
