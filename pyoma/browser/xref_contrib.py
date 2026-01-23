import collections
import io
import multiprocessing
import itertools
import re
import numpy
import pandas
import tables
import logging
import os
from pathlib import Path
from typing import List, Tuple, Set, Union
from .hogprofile.build import Pipeline, SourceProcess, BaseProfileBuilderProcess, Stage

logger = logging.getLogger(__name__)


class GeneEntries:
    def __init__(self, enrs, main):
        self.main = main
        self.enrs = enrs

    def entrynr_slices(self):
        if self.enrs[-1] - self.enrs[0] == len(self.enrs) - 1:
            return (slice(self.enrs[0], self.enrs[-1] + 1, 1),)
        slices, i0 = [], 0
        for i in range(len(self.enrs)):
            if self.enrs[i] - self.enrs[i0] == i - i0:
                continue
            slices.append(slice(self.enrs[i0], self.enrs[i - 1] + 1, 1))
            i0 = i
        slices.append(slice(self.enrs[i0], self.enrs[i] + 1, 1))
        return tuple(slices)


class SpliceVariantHelper:
    def __init__(self, h5, alt_array=None):
        if alt_array is not None:
            lookup = numpy.copy(alt_array)
            lookup[numpy.where(lookup > 0)] -= 1
        else:
            pe = h5.get_node("/Protein/Entries")
            lookup = numpy.zeros(len(pe), dtype="i4")
            for row in pe.where("AltSpliceVariant > 0"):
                lookup[row.nrow] = row["AltSpliceVariant"] - 1
        self.lookup = lookup
        self.n = len(lookup)

    def iter_genes(self):
        taken = numpy.zeros(self.n, dtype=bool)
        for i in range(self.n):
            if taken[i]:
                continue
            if self.lookup[i] == 0:
                taken[i] = True
                yield GeneEntries(numpy.array([i + 1], dtype="i4"), i + 1)
            else:
                main = self.lookup[i]
                l = self.lookup[i : i + 1000]
                res = numpy.where(l == main)[0] + i
                taken[res] = True
                yield GeneEntries(res + 1, main + 1)


class KeywordIndexer:
    def __init__(self, stopwords: Union[os.PathLike, Set[str]]):
        self.kw = collections.defaultdict(set)
        if isinstance(stopwords, set):
            self.stopwords = set(stopwords)
        else:
            self.stopwords = set()
            with open(stopwords, "rt") as fh:
                for line in fh:
                    self.stopwords.union(line.split())

    def process(self, desc, enr):
        pass


class XRefIndexHandler(BaseProfileBuilderProcess):
    def __init__(self, outfile, **kwargs):
        super().__init__(**kwargs)
        if outfile is None:
            outfile = Path(os.getenv("TMPDIR", "/tmp")) / "tmp_index.h5"
        self.outfile = outfile
        self.tmp_h5 = outfile + ".tmp"
        self.kwi = None

    def setup(self):
        self.xref_h5 = tables.open_file(self.tmp_h5, "w", filters=tables.Filters(7, complib="blosc2"))
        xref_dtype = numpy.dtype([("XRefId", "S50"), ("EntryNr", "i4"), ("XRefRow", "i4")])
        grp = self.xref_h5.create_group("/", "XRefIndex", title="auxiliary lookup tables with deduplicated xrefs")
        self.xref_idx = self.xref_h5.create_table(
            grp,
            "XRefs",
            description=xref_dtype,
            chunkshape=(16384,),
        )
        self.genenames = collections.defaultdict(list)
        self.spids = collections.defaultdict(list)
        self._buffer = []
        self.kwi = KeywordIndexer(stopwords=set([]))

    def _sort_and_store_xrefs(self):
        data = self.xref_idx.read()
        data.sort(order=["XRefId", "EntryNr"])
        self.xref_h5.close()
        self.xref_h5 = tables.open_file(self.outfile, "w", filters=tables.Filters(7, complib="blosc2"))
        xref_idx = self.xref_h5.create_table(
            "/XRefIndex", "XRefs", obj=data, expectedrows=len(data), createparents=True
        )
        xref_idx.colinstances["XRefId"].create_csindex()
        self.xref_idx = xref_idx

    def _store_names(self, dic, node, skip_short=0):
        keys = sorted(dic.keys())
        max_key_len = max(len(z) for z in keys)
        dtyp_key = numpy.dtype([("XRefId", "S{}".format(max_key_len)), ("Offset", "i4"), ("Length", "i4")])
        dtyp_lookup = numpy.dtype([("EntryNr", "i4"), ("XRefRow", "i4")])
        tab_key = self.xref_h5.create_table("/XRefIndex", node, description=dtyp_key, expectedrows=len(keys))
        tab_lookup = self.xref_h5.create_table(
            "/XRefIndex",
            node + "_lookup",
            description=dtyp_lookup,
            expectedrows=sum(len(z) for z in dic.values()),
        )
        off, skipped = 0, 0
        for k in keys:
            if len(dic[k]) < skip_short:
                logger.debug("skipping %s from index", k)
                skipped += 1
                continue
            tab_lookup.append(sorted(dic[k]))
            tab_key.append([(k, off, len(dic[k]))])
            off += len(dic[k])
        tab_lookup.flush()
        tab_key.flush()
        if len(tab_lookup) != off:
            raise Exception("Lookup table length does not match key table offsets")
        tab_key.colinstances["XRefId"].create_csindex()
        logger.info("stored {} index with {} elements. skipped {} keys".format(node, len(tab_key), skipped))

    def _add_to_buffer(self, e):
        self._buffer.append(e)
        if len(self._buffer) >= 2 * self.xref_idx.chunkshape[0]:
            self._flush()

    def _flush(self):
        self.xref_idx.append(self._buffer)
        self._buffer = []

    def add_xref(self, xref, enr, xref_row):
        xref = xref.lower()
        self._add_to_buffer((xref, enr, xref_row))

    def add_swissprot(self, id, enr, xref_row):
        self.spids[id.lower()].append((enr, xref_row))
        self.add_xref(id, enr, xref_row)

    def add_gene_name(self, id, enr, xref_row):
        self.genenames[id.lower()].append((enr, xref_row))
        self.add_xref(id, enr, xref_row)

    def add_keyword(self, enr: int, desc: str):
        self.kwi.process(desc, enr)

    def handle_input(self, item: Tuple[pandas.DataFrame, List]):
        df, desc = item
        for row in df.to_records(index=False):
            if row["XRefSource"] == 0:
                self.add_swissprot(row["XRefId"], row["EntryNr"], row["xref_row"])
            elif row["XRefSource"] in (110, 115):
                self.add_gene_name(row["XRefId"], row["EntryNr"], row["xref_row"])
            else:
                self.add_xref(row["XRefId"], row["EntryNr"], row["xref_row"])

    def finalize(self):
        self._flush()
        self._sort_and_store_xrefs()
        self._store_names(self.spids, "SwissProt")
        self._store_names(self.genenames, "GeneNames", 3)
        self.xref_h5.close()
        os.remove(self.tmp_h5)


def rem_vers(x: bytes) -> bytes:
    i = x.rfind(b".")
    if i != -1:
        suf = x[i + 1 :]
        if 1 <= len(suf) <= 2 and suf.isdigit():
            return x[:i]
    return x


class GeneGenerator(SourceProcess):
    def __init__(self, h5_path, **kwargs):
        super().__init__(**kwargs)
        self.h5_path = h5_path
        self.h5 = None
        self.splice_helper = None

    def setup(self):
        self.h5 = tables.open_file(self.h5_path)
        self.splice_helper = SpliceVariantHelper(self.h5)

    def generate_data(self):
        for gene in self.splice_helper.iter_genes():
            yield gene

    def finalize(self):
        self.h5.close()


class XRefReducer(BaseProfileBuilderProcess):
    def __init__(self, xref_path, db_path=None, **kwargs):
        super().__init__(**kwargs)
        self.xref_path = xref_path
        self.db_path = db_path

        self.xref = None
        self.db = None
        self.xref_tab = None

        # cached columns for fast access
        self._xref_eof = None
        self._prot_entries = None
        self._desc_buf = None

    def setup(self):
        self.xref = tables.open_file(self.xref_path)
        if self.db_path is not None and self.db_path != self.xref_path:
            self.db = tables.open_file(self.db_path)
        else:
            self.db = self.xref

        self.xref_tab = self.xref.get_node("/XRef")
        try:
            self._xref_eof = self.xref.get_node("/XRef_EntryNr_offset").read()
        except tables.NoSuchNodeError:
            # build xref_eof from xref_tab
            self._xref_eof = self._build_xref_eof()
        self._prot_entries = self.db.get_node("/Protein/Entries")
        self._desc_buf = self.db.get_node("/Protein/DescriptionBuffer")

    def finalize(self):
        # release cached columns
        self._xref_eof = None
        self._prot_entries = None
        self._desc_buf = None
        # close files
        if self.db != self.xref:
            self.db.close()
        self.xref.close()

    def _build_xref_eof(self):
        enr_col = self.xref_tab.col("EntryNr")
        if not (enr_col[:-1] <= enr_col[1:]).all():
            raise RuntimeError("EntryNr column is not sorted in /XRef")
        enrs = numpy.arange(enr_col[-1] + 2)
        idx = numpy.searchsorted(enr_col, enrs).astype("i4")
        return idx

    def _load_xrefs(self, slices: List[slice]) -> Tuple[numpy.ndarray | None, numpy.ndarray | None]:
        blocks = []
        row_blocks = []
        for s in slices:
            lo = self._xref_eof[s.start]
            hi = self._xref_eof[s.stop]
            if lo == hi:
                continue

            blocks.append(self.xref_tab[lo:hi])
            row_blocks.append(numpy.arange(lo, hi, dtype=numpy.int64))

        if not blocks:
            return None, None
        data = numpy.concatenate(blocks)
        rows = numpy.concatenate(row_blocks)
        return data, rows

    def _load_descriptions(self, slices: List[slice]):
        descriptions = []
        for s in slices:
            entries = self._prot_entries[s.start - 1 : s.stop - 1]
            for r in entries:
                off = r["DescriptionOffset"]
                length = r["DescriptionLength"]
                desc = self._desc_buf[off : off + length].tobytes().decode()
                descriptions.append((r["EntryNr"], desc))
        return descriptions

    def handle_input(self, gene):
        slices = gene.entrynr_slices()
        xref_data, xref_rows = self._load_xrefs(slices)
        descs = self._load_descriptions(slices)

        if xref_data is None:
            return pandas.DataFrame(), descs

        # compute is_main mask
        is_main = xref_data["EntryNr"] == gene.main

        # normalize XRefIds (lowercase, remove version suffix)
        norm = [rem_vers(x.lower()) for x in xref_data["XRefId"]]
        maxlen = max(map(len, norm))
        xref_id = numpy.array(norm, dtype=f"S{maxlen}")

        # --------------------------------------------------- sort keys
        # sort priority:
        #   Verification ASC
        #   is_main DESC
        #   XRefSource ASC
        order = numpy.lexsort(
            (
                xref_data["XRefSource"],
                ~is_main,
                xref_data["Verification"],
            )
        )

        data = xref_data[order]
        rows = xref_rows[order]
        xref_id = xref_id[order]
        is_main = is_main[order]

        # ---------------------------------------- deduplicate by XRefId
        _, idx = numpy.unique(xref_id, return_index=True)
        data = data[idx]
        rows = rows[idx]
        xref_id = xref_id[idx]
        is_main = is_main[idx]

        # ---------------------------------------- build pandas dataframe
        df = pandas.DataFrame.from_records(data)
        df["xref_row"] = rows
        df["XRefId"] = xref_id
        df["is_main"] = is_main

        return df.reset_index(drop=True), descs


def reduce_xrefs(h5_path, xref_path=None, outpath=None, nr_procs=None):
    pipeline = Pipeline()
    if nr_procs is None:
        nr_procs = multiprocessing.cpu_count()

    if xref_path is None:
        xref_path = h5_path
    pipeline.add_stage(Stage(GeneGenerator, nr_procs=1, h5_path=h5_path))
    pipeline.add_stage(Stage(XRefReducer, nr_procs=nr_procs, xref_path=xref_path, db_path=h5_path))
    pipeline.add_stage(Stage(XRefIndexHandler, nr_procs=1, outfile=outpath))
    print("setup pipeline, about to start it")
    pipeline.run()
    print("finished computing a reduced set of xrefs")
