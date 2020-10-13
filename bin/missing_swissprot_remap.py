import pyoma.browser.db
from pyoma.hpc import detect_hpc_jobarray
import gzip
import Bio.SeqIO
import logging
import os
import csv

logger = logging.getLogger(__name__)


def get_missing_recs(db, fn, format):
    all_ids = set(db.get_hdf5_handle().root.XRef.col("XRefId"))
    all_taxids = set(db.get_hdf5_handle().root.Genome.col("NCBITaxonId"))
    with gzip.open(fn, "rt") as fh:
        missing = []
        for rec in Bio.SeqIO.parse(fh, format=format):
            if int(rec.annotations["ncbi_taxid"][0]) in all_taxids:
                if rec.id.encode("utf-8") not in all_ids:
                    missing.append(rec)
    return missing


def identity_w_gaps_ignored(al1, al2):
    n = len(al1)
    match = sum(al1[k] == al2[k] or al1[k] == "_" or al1[k] == "_" for k in range(n))
    return match / n


def search_rec(db, rec):
    g = db.id_mapper["OMA"].genome_from_taxid(int(rec.annotations["ncbi_taxid"][0]))
    rng = g["EntryOff"] + 1, g["EntryOff"] + g["TotEntries"]
    match = db.seq_search.approx_search(
        str(rec.seq), compute_distance=True, entrynr_range=rng
    )
    if match is None:
        return None

    best, stats = match[0]
    if (
        stats["score"] > 1000
        and identity_w_gaps_ignored(stats["alignment"][0][0], stats["alignment"][1][0])
        > 0.8
    ):
        logger.info(
            "map {}({}) to {}: score: {}, ident: {}, kmer: {}".format(
                rec.name,
                rec.id,
                best,
                stats["score"],
                identity_w_gaps_ignored(
                    stats["alignment"][0][0], stats["alignment"][1][0]
                ),
                stats["kmer_coverage"],
            )
        )
        return best


def map_missing(db, missing, nr_procs=1):
    tab_ext = []
    src_enum, verif_enum = (
        db.id_mapper["XRef"].xref_tab.get_enum(k)
        for k in ("XRefSource", "Verification")
    )
    pInf = detect_hpc_jobarray(nr_procs)
    for rec in missing:
        if not pInf.is_my_job(rec.id):
            continue
        best = search_rec(db, rec)
        if best is not None:
            tab_ext.extend(
                [
                    (
                        best,
                        src_enum["UniProtKB/TrEMBL"],
                        rec.id.encode("utf-8"),
                        verif_enum["modified"],
                    ),
                    (
                        best,
                        src_enum["UniProtKB/SwissProt"],
                        rec.name.encode("utf-8"),
                        verif_enum["modified"],
                    ),
                ]
            )
    write_matches(pInf.modify_filename("extra-xrefs.tsv"), tab_ext)
    return tab_ext


def write_matches(fn, matches):
    with open(fn, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerows(matches)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="mapping a swissprot data file")
    parser.add_argument("--db", required=True, help="hdf5-db path")
    parser.add_argument(
        "--seqs", required=True, help="path to a file with sequences to be mapped"
    )
    parser.add_argument("--format", default="swiss", help="format of seqs file")
    parser.add_argument("--phase", choices=("filter", "map"), required=True)
    parser.add_argument("--nr-procs", "-n", default=1, type=int)
    conf = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    nr_procs = conf.nr_procs
    if nr_procs is None:
        nr_procs = int(os.getenv("NR_PROCESSES", "1"))
    db = pyoma.browser.db.Database(conf.db)
    if conf.phase == "filter":
        missing = get_missing_recs(db, conf.seqs, conf.format)
        with open(conf.seqs + ".missings", "wt") as fh:
            Bio.SeqIO.write(missing, fh, format=conf.format)
    elif conf.phase == "map":
        with open(conf.seqs + ".missings", "rt") as fh:
            missing = list(Bio.SeqIO.parse(fh, format=conf.format))
        tab_ext = map_missing(db, missing, nr_procs=nr_procs)
    db.close()
