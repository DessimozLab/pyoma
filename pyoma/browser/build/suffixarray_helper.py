import logging
import tables
import numpy
from tqdm import tqdm

logger = logging.getLogger(__name__)


def build_filtered_sa(
    sa_h5: tables.CArray,  # permanent suffix array (CArray in /Protein/SequenceIndex)
    delimiters_sorted: numpy.ndarray,  # sorted array of entry delimiters (first nr_entries of SA)
    n_seq: int,  # length of concatenated sequence buffer
    nr_entries: int,  # number of sequences
    k: int,  # k-mer size for filtering
    dtype_sa: numpy.dtype,  # dtype of suffix array positions (from sa_h5.dtype)
    dtype_enr: numpy.dtype,  # dtype for entry numbers (uint32 or uint64)
    h5_out: tables.File,  # output HDF5 file handle (temporary or permanent)
):
    """
    Stream over the raw suffix array (sa_h5) and build filtered SA, original SA positions,
    and entry indices in SA-order, all stored in h5_out.

    Returns:
        (sa_f, sa_origpos, idx_saorder) as PyTables EArrays
    """
    # add sentinel delimiter at the end (ensure every position has a "next delimiter")
    if delimiters_sorted[-1] < n_seq - 1:
        delimiters_sorted = numpy.concatenate(
            [delimiters_sorted, numpy.array([n_seq - 1], dtype=delimiters_sorted.dtype)]
        )

    sa_f = h5_out.create_earray(
        "/",
        "SequenceIndex",
        atom=tables.Atom.from_dtype(dtype_sa),
        shape=(0,),
        expectedrows=n_seq,
        title=f"Filtered SA (k={k})",
    )
    sa_origpos = h5_out.create_earray(
        "/",
        "SequenceIndexOrigPos",
        atom=tables.Atom.from_dtype(dtype_sa),
        shape=(0,),
        expectedrows=n_seq,
        title=f"Original SA indices (k={k})",
    )
    idx_saorder = h5_out.create_earray(
        "/",
        "IndexInSAOrder",
        atom=tables.Atom.from_dtype(numpy.dtype(dtype_enr)),
        shape=(0,),
        expectedrows=n_seq,
        title=f"Entry numbers aligned to filtered SA order (k={k})",
    )

    logger.info(f"Streaming SA to filter with k={k}...")
    sa_chunk_read = 10_000_000
    for start in tqdm(range(0, len(sa_h5), sa_chunk_read), desc=f"Filter SA (k={k})"):
        stop = min(len(sa_h5), start + sa_chunk_read)
        chunk_sa = sa_h5[start:stop]

        # Vectorized lookup of nearest delimiter >= p
        idxs = numpy.searchsorted(delimiters_sorted, chunk_sa, side="left")  # int[chunk]
        e = delimiters_sorted[idxs]  # uint64[chunk]

        # Keep iff (e - p) >= k
        keep = (e - chunk_sa) >= k
        if not numpy.any(keep):
            continue

        kept_sa = chunk_sa[keep]
        kept_origpos = (start + numpy.nonzero(keep)[0]).astype(dtype_sa)
        kept_entries = (idxs[keep] + 1).astype(dtype_enr)  # 1-based

        sa_f.append(kept_sa)
        sa_origpos.append(kept_origpos)
        idx_saorder.append(kept_entries)

    # flush to disk
    sa_f.flush()
    sa_origpos.flush()
    idx_saorder.flush()

    return sa_f, sa_origpos, idx_saorder


def kmer_codes_for_positions(P, k, seqs_np, dtype_sa, map256, alphabet_size):
    """
    Vectorized equivalent of [KmerEncoder.decode(seqs[p:p+k]) for p in P].

    P: np.ndarray of start positions
    k: int,   k-mer size
    seqs_np:  np.ndarray view of the sequence buffer
    dtype_sa: np.dtype for the position data
    map256:   map from character to uint8
    alphabet_size: int, size of the alphabet (21 for AA, 5 for DNA)
    Returns: codes (dtype_sa), valid_mask (bool)
    """
    P = P.astype(dtype_sa, copy=False)
    codes = numpy.zeros(len(P), dtype=dtype_sa)
    good = numpy.ones(len(P), dtype=bool)

    for t in range(k):
        b = seqs_np[P + t]
        v = map256[b]
        bad = v == 255
        good &= ~bad
        codes = codes * alphabet_size + v.astype(dtype_sa)

    codes[~good] = numpy.iinfo(dtype_sa).max  # sentinel
    return codes, good
