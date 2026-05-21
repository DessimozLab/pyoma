from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import tables

from .exceptions import DBConsistencyError

if TYPE_CHECKING:
    from .db import Database


@dataclass
class StructureInfo:
    seq_3di: bytes
    sequence: bytes  # AA sequence, same length as seq_3di
    source: str  # "AlphaFold" | "ProstT5"


class StructureDB:
    def __init__(self, db: Database, path: os.PathLike):
        self.h5 = tables.open_file(path, "r")
        self._db = db
        self._index: tables.Table = self.h5.get_node("/index")
        self._seq_3di: tables.EArray = self.h5.get_node("/sequences_3di")
        self._source_enum = self._index.get_enum("Source")
        db.register_on_close(self.close)

    def close(self):
        self.h5.close()
        self._db.unregister_on_close(self.close)

    def get(self, entry) -> Optional[StructureInfo]:
        entry_nr = int(entry["EntryNr"])
        row = self._index[entry_nr - 1]
        if row["EntryNr"] != entry_nr:
            rows = self._index.read_where("EntryNr == entry_nr", condvars={"entry_nr": entry_nr})
            if len(rows) == 0:
                return None
            if len(rows) > 1:
                raise DBConsistencyError(f"Expected exactly one row for EntryNr {entry_nr}, but found {len(rows)}")
            row = rows[0]
        if row["Length_3DI"] == 0:
            return None

        off, length = int(row["Offset_3DI"]), int(row["Length_3DI"])
        seq_3di = self._seq_3di[off : off + length - 1].tobytes()
        source = self._source_enum(int(row["Source"]))
        sequence = self._db.get_sequence(entry)
        return StructureInfo(seq_3di=seq_3di, sequence=sequence, source=source)
