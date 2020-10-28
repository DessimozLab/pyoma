#!/usr/bin/env python
import csv
import pandas
import pyoma.browser.db


def go_annotations(db, species=None, evidence=None, reference=None, idtype=None):
    if not isinstance(db, pyoma.browser.db.Database):
        raise TypeError("db must be a hdf5 oma database")
    if species is not None:
        start, stop = db.id_mapper["OMA"].genome_range(species)
    else:
        start, stop = (1, db.id_resolver.max_entry_nr + 1)
    go_df = db.get_gene_ontology_annotations(start, stop=stop, as_dataframe=True)

    if evidence is not None:
        if isinstance(evidence, str):
            go_df = go_df.loc[go_df["Evidence"] == evidence]
        elif isinstance(evidence, (list, tuple)):
            go_df = go_df.loc[go_df["Evidence"].isin(evidence)]

    if reference is not None:
        if isinstance(reference, str):
            go_df = go_df.loc[go_df["DB:Reference"] == reference]
        elif isinstance(reference, (list, tuple)):
            go_df = go_df.loc[go_df["DB:Reference"].isin(reference)]

    if idtype is not None:
        xref_tab = db.db.get_node("/XRef")
        xrefs = pandas.DataFrame(
            xref_tab.read_where(
                "(EntryNr > {}) & (EntryNr <= {}) & (XRefSource == {})".format(
                    start, stop, xref_tab.get_enum("XRefSource")[idtype]
                )
            )
        )
        xrefs.drop(columns=["XRefSource", "Verification"], inplace=True)
        xrefs["XRefId"] = xrefs["XRefId"].str.decode("utf-8")
        xrefs["OmaID"] = xrefs["EntryNr"].apply(db.id_mapper["Oma"].map_entry_nr)
        new_df = go_df.merge(
            xrefs, how="left", left_on="DB_Object_ID", right_on="OmaID"
        )
        new_df.drop(columns=["DB_Object_ID"], inplace=True)
        new_df[idtype] = new_df["XRefId"]
        new_df.rename(columns={"XRefId": "DB_Object_ID"}, inplace=True)
        go_df = new_df
    else:
        go_df["OmaID"] = go_df["DB_Object_ID"]
    return go_df


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Dump GO Annotations from an hdf5 file"
    )
    parser.add_argument(
        "--genome",
        "-g",
        help="Genome id for which annotations should be returned. "
        "By default annotations for all genomes will be returned",
    )
    parser.add_argument(
        "--only-oma",
        action="store_true",
        default=False,
        help="Return only annotations produced by the OMA function annotation pipeline",
    )
    parser.add_argument(
        "--xref-type",
        help="CrossReference ID used in output. The value must be matching one of the "
        "enum values in the XRef table, e.g. SourceAC",
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
    conf.ref = "OMA_Fun:001" if conf.only_oma else None

    db = pyoma.browser.db.Database(conf.db)
    annotations = go_annotations(
        db, species=conf.genome, reference=conf.ref, idtype=conf.xref_type
    )

    cols = ["OmaID", "GO_ID", "Evidence", "DB:Reference"]
    if conf.xref_type is not None:
        cols.insert(1, conf.xref_type)

    csv_writer = csv.DictWriter(conf.out, cols, extrasaction="ignore", delimiter="\t")
    csv_writer.writerows(annotations.to_dict("records"))
