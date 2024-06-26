#!/usr/bin/env python3

import pyoma.browser.db
from pyoma.common import auto_open
import numpy
import pandas as pd


def get_assignments(db):
    df = pd.DataFrame(
        numpy.fromiter(
            ((row["EntryNr"], row["OmaHOG"].decode()) for row in db.get_hdf5_handle().root.Protein.Entries),
            dtype=[("EntryNr", "i4"), ("OmaHOG", "S255")],
        )
    )
    df["OmaID"] = df["EntryNr"].apply(db.id_mapper["Oma"].map_entry_nr)
    df["OmaHOG"] = df["OmaHOG"].apply(bytes.decode)
    assert df["EntryNr"].equals(pd.Series(data=range(1, len(df) + 1), dtype="i4"))
    return df


def load_xref(db, idtype):
    xref_tab = db.get_hdf5_handle().get_node("/XRef")
    enum = xref_tab.get_enum("XRefSource")
    sources = [enum[z] for z in idtype]
    cond = " | ".join(["(XRefSource == {})".format(z) for z in sources])
    df = pd.DataFrame(
        numpy.fromiter(
            ((row["EntryNr"], row["XRefSource"], row["XRefId"].decode()) for row in xref_tab.where(cond)),
            dtype=[("EntryNr", "i4"), ("Source", "i1"), ("XRefId", "S100")],
        )
    )
    df["Source"] = df["Source"].apply(lambda s: enum(s))
    df["XRefId"] = df["XRefId"].apply(bytes.decode)
    wide = df.pivot_table(index="EntryNr", columns="Source", values="XRefId", aggfunc=lambda s: " _|_ ".join(s))
    wide = wide.reset_index(level="EntryNr")
    return wide


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Dump HOG ids for all protein in OMA")
    parser.add_argument("--out", required=True, help="path to output file")
    parser.add_argument("db", help="Path to input hdf5 file")
    parser.add_argument("--idtype", default=[], nargs="*", help="list of idtypes to include in the output")
    parser.add_argument("--format", choices=("tsv", "darwin"), default="tsv", help="outout format")
    conf = parser.parse_args()

    db = pyoma.browser.db.Database(conf.db)
    df = get_assignments(db)
    df = df.set_index("EntryNr")[["OmaID", "OmaHOG"]]
    if len(conf.idtype) > 0:
        xrefs = load_xref(db, conf.idtype)
        df = df.join(xrefs.set_index("EntryNr"), on="EntryNr")

    if conf.format == "tsv":
        df.to_csv(
            conf.out,
            index=False,
            compression="infer",
            sep="\t",
            # columns=["OmaID", "OmaHOG"],
        )

    if conf.format == "darwin":
        with auto_open(conf.out, "wt", encoding="utf-8") as fh:
            fh.write("HOG_IDs := [\n")
            for hogid in df["OmaHOG"]:
                fh.write("'{}',\n".format(hogid))
            fh.write("NULL]:\n")
