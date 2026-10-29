from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterator, List, Optional, Tuple

import numpy
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
        row = self._index[entry_nr - 1] if entry_nr <= self._index.nrows else None
        if row is None or row["EntryNr"] != entry_nr:
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

    def iter_seq_3di(
        self, start: int = 1, stop: Optional[int] = None, block_size: int = 50_000, max_gap: int = 1 << 20
    ) -> Iterator[Tuple[int, bytes]]:
        """Yield (entry_nr, seq_3di) ordered by EntryNr for all entries in [start, stop) with a 3Di sequence.

        Unlike :meth:`get`, the index and the 3Di buffer are read blockwise
        (``block_size`` entries at a time), which is much faster when iterating
        over many entries. Buffer ranges that are less than ``max_gap`` bytes
        apart are fetched with a single read."""
        if stop is None:
            stop = self._db.nr_proteins + 1
        for lo in range(start, stop, block_size):
            rows = self._index_rows(lo, min(lo + block_size, stop))
            rows = rows[rows["Length_3DI"] > 0]
            if len(rows) > 0:
                yield from zip(rows["EntryNr"].tolist(), self._read_seq_3di(rows, max_gap))

    def _index_rows(self, lo: int, hi: int) -> numpy.ndarray:
        """index rows of entries in [lo, hi), sorted by EntryNr"""
        rows = self._index[lo - 1 : hi - 1]
        if len(rows) == hi - lo and numpy.array_equal(rows["EntryNr"], numpy.arange(lo, hi)):
            return rows
        rows = self._index.read_where("(EntryNr >= lo) & (EntryNr < hi)", condvars={"lo": lo, "hi": hi})
        rows.sort(order="EntryNr")
        if len(rows) > 1 and (numpy.diff(rows["EntryNr"].astype(numpy.int64)) == 0).any():
            raise DBConsistencyError(f"Expected at most one row per EntryNr in [{lo}, {hi})")
        return rows

    def _read_seq_3di(self, rows: numpy.ndarray, max_gap: int) -> List[bytes]:
        offsets = rows["Offset_3DI"].astype(numpy.int64)
        lengths = rows["Length_3DI"].astype(numpy.int64)
        order = numpy.argsort(offsets, kind="stable")
        seqs = [b""] * len(rows)
        i = 0
        while i < len(order):
            # extend the segment while the next sequence starts close enough to its end
            seg_start, seg_end, j = offsets[order[i]], offsets[order[i]] + lengths[order[i]], i + 1
            while j < len(order) and offsets[order[j]] - seg_end <= max_gap:
                seg_end = max(seg_end, offsets[order[j]] + lengths[order[j]])
                j += 1
            buf = self._seq_3di[seg_start:seg_end].tobytes()
            for k in order[i:j]:
                o = offsets[k] - seg_start
                seqs[k] = buf[o : o + lengths[k] - 1]  # stored length includes the terminator
            i = j
        return seqs
