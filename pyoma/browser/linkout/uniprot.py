"""Generate crosslinks for UniProtKB pointing to OMA."""

import csv
import logging

from ..db import Database
from pyoma.common import auto_open

logger = logging.getLogger(__name__)


class OmaGroupMapper(object):
    def __init__(self, db):
        grp2finger = {}
        for row in db.get_hdf5_handle().get_node("/OmaGroups/MetaData"):
            if row["Fingerprint"] != b"n/a":
                grp2finger[int(row["GroupNr"])] = row["Fingerprint"].decode()
        self.grp2finger = grp2finger
        self.db = db

    def map(self, entry_nrs):
        for enr, xref in entry_nrs:
            e = self.db.entry_by_entry_nr(enr)
            try:
                fingerprint = self.grp2finger[e["OmaGroup"]]
                yield xref, fingerprint
            except KeyError:
                pass


def iter_uniprot_entries(xref_tab):
    enum = xref_tab.get_enum("XRefSource")
    tr = enum["UniProtKB/TrEMBL"]
    for row in xref_tab:
        if row["XRefSource"] == tr:
            yield int(row["EntryNr"]), row["XRefId"].decode()


def iter_uniprot_omagrp_mapping(hdf5):
    with Database(hdf5) as db:
        grp_mapper = OmaGroupMapper(db)
        up_entries = iter_uniprot_entries(db.get_hdf5_handle().get_node("/XRef"))
        yield from grp_mapper.map(up_entries)


def write_omagrp_xref_map(fh, hdf5):
    csv_writer = csv.writer(fh, dialect="excel-tab")
    csv_writer.writerows(iter_uniprot_omagrp_mapping(hdf5))
