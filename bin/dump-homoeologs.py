#!/usr/bin/env python
import csv
import pandas
import tables
from pyoma.browser.homoeologs import HomeologsConfidenceCalculator


class HomoeologExtractor(HomeologsConfidenceCalculator):
    fields_to_load = [
        "EntryNr1",
        "EntryNr2",
        "Confidence",
        "Distance",
        "SyntenyConservationLocal",
    ]

    def add_more_protein_features(self, df):
        entry_tab = self.h5_handle.root.Protein.Entries
        entries = pandas.DataFrame(
            entry_tab.read_where("(EntryNr >= {}) & (EntryNr <= {})".format(self.genome_range[0], self.genome_range[1]))
        )
        entries = entries[["EntryNr", "Chromosome", "SubGenome"]]
        for c in ("Chromosome", "SubGenome"):
            entries[c] = entries[c].str.decode("utf-8")
        entries = entries.assign(
            OmaID=entries["EntryNr"].apply(lambda x: "{:s}{:05d}".format(self.genome, x - self.genome_range[0] + 1))
        )

        new_df = pandas.merge(
            pandas.merge(df, entries, how="left", left_on="EntryNr1", right_on="EntryNr"),
            entries,
            how="left",
            left_on="EntryNr2",
            right_on="EntryNr",
        )
        new_df.drop(columns=["EntryNr_x", "EntryNr_y"], inplace=True, errors="ignore")
        new_df.rename(
            columns={
                "OmaID_x": "OmaID1",
                "OmaID_y": "OmaID2",
                "Chromosome_x": "Chromosome1",
                "Chromosome_y": "Chromosome2",
                "SubGenome_x": "SubGenome1",
                "SubGenome_y": "SubGenome2",
            },
            inplace=True,
        )
        return new_df


def homoeologs_of_genome(db, genome):
    extractor = HomoeologExtractor(db, genome)
    data = extractor.augment_ids(extractor.relations_df, keep="agg")
    data = extractor.add_more_protein_features(data)
    return data


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Dump Homoeolog Information from an hdf5 file")
    parser.add_argument(
        "--genome",
        "-g",
        help="Genome id for which annotations should be returned. "
        "By default annotations for all genomes will be returned",
    )
    parser.add_argument(
        "--out",
        "-o",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file path, defaults to stdout",
    )
    parser.add_argument("db", help="Path to the hdf5 database file")
    conf = parser.parse_args()

    cols = [
        "OmaID1",
        "SourceID1",
        "Chromosome1",
        "SubGenome1",
        "OmaID2",
        "SourceID2",
        "Chromosome2",
        "SubGenome2",
        "Confidence",
        "Distance",
        "SyntenyConservationLocal",
    ]
    csv_writer = csv.DictWriter(conf.out, cols, extrasaction="ignore", delimiter="\t")
    csv_writer.writeheader()

    if conf.genome is None:
        with tables.open_file(conf.db) as h5:
            genomes = [x["UniProtSpeciesCode"].decode() for x in h5.root.Genome.read_where("IsPolyploid==True")]
    else:
        genomes = [conf.genome]

    for g in genomes:
        data = homoeologs_of_genome(conf.db, g)
        csv_writer.writerows(data.to_dict("records"))
