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
        self.xref_h5 = tables.open_file(self.tmp_h5, "w", filters=tables.Filters(6, complib="blosc"))
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
        self.kwi = KeywordIndexer()

    def _sort_and_store_xrefs(self):
        data = self.xref_idx.read()
        data.sort(order=["XRefId", "EntryNr"])
        self.xref_h5.close()
        self.xref_h5 = tables.open_file(self.outfile, "w", filters=tables.Filters(6, complib="blosc"))
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


RE_VERS = re.compile(rb"(?P<base>[a-z0-9_.-]+)\.\d{1,2}$")


def rem_vers(x):
    m = RE_VERS.match(x)
    if m is not None:
        return m.group("base")
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
    def __init__(self, h5_path, **kwargs):
        super().__init__(**kwargs)
        self.h5_path = h5_path
        self.h5 = None
        self.xref_tab = None
        self.xref_eof = None

    def setup(self):
        self.h5 = tables.open_file(self.h5_path)
        self.xref_tab = self.h5.get_node("/XRef")
        try:
            self.xref_eof = self.h5.get_node("/XRef_EntryNr_offset").read()
        except tables.NoSuchNodeError:
            pass

    def finalize(self):
        self.h5.close()

    def _load_xrefs_with_entry_offset(self, gene):
        xrefs = pandas.DataFrame(
            numpy.hstack(
                list(
                    map(
                        lambda s: self.xref_tab[self.xref_eof[s.start] : self.xref_eof[s.stop]],
                        gene.entrynr_slices(),
                    )
                )
            )
        )
        xrefs["xref_row"] = pandas.Series(
            itertools.chain.from_iterable(
                map(
                    lambda s: range(self.xref_eof[s.start], self.xref_eof[s.stop]),
                    gene.entrynr_slices(),
                )
            )
        )
        return xrefs

    def _load_xrefs_with_where_cond(self, gene):
        query = " | ".join(
            ["((EntryNr >= {}) & (EntryNr < {}))".format(s.start, s.stop) for s in gene.entrynr_slices()]
        )
        dtype = [("xref_row", "i8")] + self.xref_tab.dtype.descr
        buf = io.BytesIO()
        for row in self.xref_tab.where(query):
            buf.write(row.nrow.tobytes())
            buf.write(row.fetch_all_fields().tobytes())
        data = numpy.frombuffer(buf.getbuffer(), dtype=dtype)
        return pandas.DataFrame(data)

    def _load_descriptions(self, gene):
        query = " | ".join([f"((EntryNr >= {s.start}) & (EntryNr < {s.stop}))" for s in gene.entrynr_slices()])
        descriptions = []
        desc_buf = self.h5.get_node("/Protein/DescriptionBuffer")
        for row in self.h5.get_node("/Protein/Entries").where(query):
            desc = (
                desc_buf[row["DescriptionOffset"] : row["DescriptionOffset"] + row["DescriptionLength"]]
                .tobytes()
                .decode()
            )
            descriptions.append((row["EntryNr"], desc))
        return descriptions

    def handle_input(self, gene):
        if self.xref_eof is None:
            xrefs = self._load_xrefs_with_where_cond(gene)
        else:
            xrefs = self._load_xrefs_with_entry_offset(gene)
        xrefs["is_main"] = xrefs["EntryNr"] == gene.main

        # transformations
        xrefs["XRefId"] = xrefs["XRefId"].apply(bytes.lower).apply(rem_vers)

        # sort such that first element per xrefid is the one we want to keep in the index
        sorted_xrefs = xrefs.sort_values(["Verification", "is_main", "XRefSource"], ascending=[True, False, True])
        res = sorted_xrefs.groupby("XRefId").first().reset_index()
        desc = self._load_descriptions(gene)
        return res, desc


def reduce_xrefs(h5_path, outpath=None, nr_procs=None):
    pipeline = Pipeline()
    if nr_procs is None:
        nr_procs = multiprocessing.cpu_count()

    pipeline.add_stage(Stage(GeneGenerator, nr_procs=1, h5_path=h5_path))
    pipeline.add_stage(Stage(XRefReducer, nr_procs=nr_procs, h5_path=h5_path))
    pipeline.add_stage(Stage(XRefIndexHandler, nr_procs=1, outfile=outpath))
    print("setup pipeline, about to start it")
    pipeline.run()
    print("finished computing a reduced set of xrefs")
