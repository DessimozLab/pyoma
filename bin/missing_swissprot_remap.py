import pyoma.browser.db
import gzip
import Bio.SeqIO
import logging
import concurrent.futures

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


def open_db(filename):
    global h5db
    h5db = pyoma.browser.db.Database(filename)


def search_rec(rec):
    g = h5db.id_mapper["OMA"].genome_from_taxid(int(rec.annotations["ncbi_taxid"][0]))
    rng = g["EntryOff"] + 1, g["EntryOff"] + g["TotEntries"]
    match = h5db.seq_search.approx_search(
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


def map_missing(db, missing):
    tab_ext = []
    src_enum, verif_enum = (
        db.id_mapper["XRef"].xref_tab.get_enum(k)
        for k in ("XRefSource", "Verification")
    )
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=32, initializer=open_db, initargs=(db.db.filename,)
    ) as executor:
        for rec, best in zip(missing, executor.map(search_rec, missing)):
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
    return tab_ext


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="mapping a swissprot data file")
    parser.add_argument("--db", required=True, help="hdf5-db path")
    parser.add_argument(
        "--seqs", required=True, help="path to a file with sequences to be mapped"
    )
    parser.add_argument("--format", default="swiss", help="format of seqs file")
    conf = parser.parse_args()

    db = pyoma.browser.db.Database(conf.db)
    missing = get_missing_recs(db, conf.seqs, conf.format)
    tab_ext = map_missing(db, missing)
    db.close()

    import csv

    with open("extra_xrefs.tsv", "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerows(tab_ext)
