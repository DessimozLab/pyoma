from typing import Optional, Dict
import logging
import os
import tempfile

import numpy
import tables
from tqdm import tqdm

from .suffixarray_helper import build_filtered_sa, kmer_codes_for_positions
from ..db import Database, SequenceSearch, read_table_where
from ..KmerEncoder import KmerEncoder, DIGITS_AA

logger = logging.getLogger(__name__)


def find_fingerprints_streaming(
    db_path: str,
    seqs_path: Optional[str] = None,
    suffix_path: Optional[str] = None,
    k: int = 7,
    batch_size: int = 1_000_000,
) -> Dict[int, str]:
    """
    Find one fingerprint (kmer length k) per OMA group using an streaming approach.

    Returns: dict {oma_group_id: kmer_str} (groups without fingerprints omitted)
    """
    # --- Preparations -------------------------------------------------------
    kmers = KmerEncoder(k, is_protein=True)
    tot_kmers = len(kmers)

    with tables.open_file(db_path, mode="r") as h5_db:
        # Build entry_nr -> group_id mapping (numpy array indexed by EntryNr)
        pe_tab: tables.Table = h5_db.get_node("/Protein/Entries")
        # Find max entry number to size the array (entries are 1-based)
        # Build mapping using vectorized fill
        nr_entries = int(pe_tab[-1]["EntryNr"])

        # for k<=7 on protein 21^7 ~1.8e9
        dtype_kmer = numpy.dtype(numpy.uint32 if (tot_kmers < numpy.iinfo(numpy.uint32).max) else numpy.uint64)

        # Fill mapping
        entry_to_group = numpy.full(nr_entries + 1, 0, dtype=numpy.int32)  # default -1
        for row in pe_tab:
            entry_to_group[row["EntryNr"]] = row["OmaGroup"]
            if row["AltSpliceVariant"] > 0 and row["AltSpliceVariant"] != row["EntryNr"]:
                entry_to_group[row["EntryNr"]] = -1  # ignore alt splice variants

        seqs = (
            numpy.memmap(seqs_path, dtype=numpy.uint8, mode="r")
            if seqs_path
            else h5_db.get_node("/Protein/SequenceBuffer")[:]
        )
    seqs_np = numpy.frombuffer(seqs, dtype=numpy.uint8)
    n_seq = len(seqs_np)
    # --- End preparations -------------------------------------------------------

    # -----------------------------
    # Stream-filter SA → SA_filtered, SA_origpos, and build IndexInSAOrder
    # -----------------------------
    with tables.open_file(suffix_path or db_path, mode="r") as h5_sa:
        sa_h5: tables.CArray = h5_sa.get_node("/Protein/SequenceIndex")
        dtype_sa = sa_h5.dtype
        delimiters_sorted = sa_h5[:nr_entries]

        # Build filtered SA for this k
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_tmp = tables.open_file(
                os.path.join(tmpdir, "sa_kmerfilter.h5"),
                mode="w",
                filters=tables.Filters(complevel=3, complib="blosc2", bitshuffle=True),
            )
            sa_f, sa_origpos, idx_saorder = build_filtered_sa(
                sa_h5, delimiters_sorted, n_seq, nr_entries, k, dtype_sa, numpy.uint32, h5_tmp
            )

            # Prepare alphabet map (AA)
            map256 = numpy.full(256, 255, dtype=numpy.uint8)
            for i, aa in enumerate(DIGITS_AA):
                map256[ord(aa)] = i
            alphabet_size = len(DIGITS_AA)

            logger.info(f"Streaming filtered SA for k={k} to find fingerprints...")
            fingerprints = {}
            L, ii = len(sa_f), 0
            # carry-over for cross-batch runs
            prev_code, prev_groups = None, None
            t = tqdm(total=max(L, 0), desc="K-mer fingerprints")
            while ii < L:
                jj = min(ii + batch_size, L)
                P = sa_f[ii:jj]
                entries = idx_saorder[ii:jj]
                codes, good = kmer_codes_for_positions(P, k, seqs_np, dtype_sa, map256, alphabet_size)

                # Find runs of equal k-mers (SA is sorted → contiguous)
                change = numpy.nonzero(codes[1:] != codes[:-1])[0] + 1
                run_starts = numpy.concatenate(([0], change))
                run_ends = numpy.concatenate((change, [codes.size]))

                for a, b in zip(run_starts, run_ends):
                    code_int = int(codes[a])
                    if code_int == numpy.iinfo(dtype_sa).max:
                        continue  # invalid k-mer (e.g., with X)
                    entry_ids = entries[a:b]
                    groups = numpy.unique(entry_to_group[entry_ids])
                    if prev_code is not None and code_int == prev_code:
                        # continuation of same k-mer from previous batch
                        groups = numpy.union1d(prev_groups, groups)
                        prev_code = None  # fully handled now
                        prev_groups = None

                    # If we're at the *last* run in this batch and it might continue,
                    # save it for the next batch to merge later.
                    is_last_run = b == len(codes)
                    if is_last_run and jj < L:
                        prev_code = code_int
                        prev_groups = groups
                        continue  # defer handling

                    # fully formed run (safe to analyze)
                    if len(groups) == 1 or (len(groups) == 2 and -1 in groups):
                        g = int(groups[groups != -1][0])
                        if g > 0 and g not in fingerprints:
                            fingerprints[g] = kmers.encode(code_int).decode()
                t.update(jj - ii)
                ii = jj

            # handle leftover carried run (if any)
            if prev_code is not None and (len(prev_groups) == 1 or (len(prev_groups) == 2 and -1 in prev_groups)):
                g = int(groups[groups != -1][0])
                if g > 0 and g not in fingerprints:
                    fingerprints[g] = kmers.encode(code_int).decode("ascii")
    for og in range(1, entry_to_group.max() + 1):
        if og not in fingerprints:
            fingerprints[og] = "n/a"
    return fingerprints
