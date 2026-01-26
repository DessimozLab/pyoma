import collections
import itertools
import os
import logging
import sys
import time
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from functools import lru_cache
from logging.handlers import QueueHandler
from multiprocessing import Queue
from multiprocessing.util import Finalize
from typing import List, Set, Optional, Iterable, Tuple
import re
import numpy
import tables

from .map import merge_sorted_h5_table
from ..builder import DBBuilder, BufferedTableWriter
from ..function_propagation_omagroup import FunctionPredictor
from ...convert import create_index_for_columns, sort_table
from ...tablefmt import GeneOntologyTable
from ...geneontology import GeneOntology, OntologyParser, AnnotationParser, GOA_Annotation, AnnotationFilter

logger = logging.getLogger(__name__)


# ------------------------------------------------------
# Generic Table Storage Manager (HDF5)
# ------------------------------------------------------
class TableStorageManager:
    """Unified HDF5 manager with buffered tables."""

    def __init__(self, filename, mode="w", buffer_size=500_000):
        self.filename = filename
        self.mode = mode
        self.buffer_size = buffer_size
        self.h5: Optional[tables.File] = None
        self.tables = {}  # table_name -> tables.Table
        self.buffers = {}  # table_name -> BufferedTableWriter

    def __enter__(self):
        self.h5 = tables.open_file(
            self.filename, mode=self.mode, filters=tables.Filters(complevel=7, complib="blosc2", fletcher32=True)
        )
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for buf in self.buffers.values():
            buf.flush()
        if self.h5:
            self.h5.flush()
            self.h5.close()
            self.h5 = None

    def create_table(self, path: str, description, sort_key=None, expectedrows=1_000_000):
        root, name = os.path.split(path)
        table = self.h5.create_table(root, name, description, expectedrows=expectedrows, createparents=True)
        self.tables[name] = table
        self.buffers[name] = BufferedTableWriter(table, sort_key=sort_key, buffer_size=self.buffer_size)
        return table

    def add_rows(self, table_name: str, rows: Iterable[Tuple]):
        self.buffers[table_name].add(rows)

    def flush_table(self, table_name: str):
        self.buffers[table_name].flush()

    def create_index(self, table_name: str, cols: List[str]):
        self.tables[table_name].flush()
        create_index_for_columns(self.tables[table_name], *cols)


def bulk_copy_table(source_tab: tables.Table, target_tab: tables.Table, chunk: int = 1_000_000):
    # 1. remove all indexes from taget_tab
    for c in target_tab.colnames:
        if target_tab.colindexed[c]:
            target_tab.colinstances[c].remove_index()

    # 2. Copy data
    for start in range(0, source_tab.nrows, chunk):
        stop = min(start + chunk, source_tab.nrows)
        target_tab.append(source_tab[start:stop])
        logger.debug(f"copied chunk {start}:{stop} of {source_tab.nrows} rows")
    target_tab.flush()


class GeneOntologyManager:
    def __init__(self, fpath, annotation_path, ontology_path):
        self.final_path = fpath
        self.annotation_path = annotation_path
        self.ontology_path = ontology_path

        self.phase = -1  # 0 -> accept ontology, 1 -> collect annotations, 2 -> infer annotations
        self.collect_file = self.final_path + ".collect"
        self.inference_file = self.final_path + ".infer"
        self.storage: Optional[TableStorageManager] = None
        self.go = None
        self._anno_tab_name = None
        self._go_buf = []
        self._inf_buf = []

    def __enter__(self):
        self.storage = TableStorageManager(self.collect_file, mode="w")
        self.storage.__enter__()
        self.storage.create_table(self.annotation_path, GeneOntologyTable, sort_key=["EntryNr", "TermNr", "Evidence"])
        self._anno_tab_name = self.annotation_path.split("/")[-1]
        self.phase = 0
        return self

    def add_ontology(self, obo_path):
        assert self.phase == 0
        # check that ontology file is not broken. if we can build it, it should be ok
        self.go = GeneOntology(OntologyParser(obo_path))
        self.go.parse()

        with open(obo_path, "rb") as fh:
            go_obo = fh.read()

        with tables.open_file(
            self.final_path, "w", filters=tables.Filters(complevel=7, complib="blosc2", fletcher32=True)
        ) as h5:
            root, name = os.path.split(self.ontology_path)
            obo = h5.create_carray(
                root,
                name,
                title="Gene ontology hierarchy definition",
                createparents=True,
                obj=numpy.ndarray(len(go_obo), buffer=go_obo, dtype=tables.StringAtom(1)),
            )
            obo.set_attr("ontology_release", self._get_obo_version(obo))
        self.phase = 1

    def flush(self):
        self.storage.flush_table(self._anno_tab_name)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.phase == 0:
            return
        if self.phase == 1:
            # we do no function prediction, but we need to copy over from collect-file
            self.switch_collect_to_inference_phase(create_index=False)
        # now, we are for sure in phase 2
        if self.phase == 2:
            self.flush()
            self.storage.create_index(self._anno_tab_name, ["EntryNr", "TermNr", "Evidence"])
            sort_table(
                self.storage.tables[self._anno_tab_name], col_order=["EntryNr", "TermNr", "Evidence", "Reference"]
            )
            self.storage.__exit__(exc_type, exc_val, exc_tb)

        # merge collect and inference files
        with TableStorageManager(self.final_path, mode="a") as storage:
            storage.create_table(self.annotation_path, GeneOntologyTable, sort_key=["EntryNr", "TermNr", "Evidence"])
            with tables.open_file(self.inference_file, "r") as h5_infer, tables.open_file(
                self.collect_file, "r"
            ) as h5_collect:
                merge_sorted_h5_table(
                    [h5_collect, h5_infer],
                    table_path=self.annotation_path,
                    sort_columns=["EntryNr", "TermNr", "Evidence", "Reference"],
                    dupl_subset_columns=["EntryNr", "TermNr", "Evidence", "Reference"],
                    storer_callback=lambda rows: storage.add_rows(self._anno_tab_name, rows),
                    enr_batch_size=1000,
                )

            storage.create_index(self._anno_tab_name, ["EntryNr", "TermNr", "Evidence"])

    def _get_obo_version(self, obo_arr):
        header = obo_arr[0:1000].tobytes()
        rel_info = re.search(rb"data-version:\s*(?P<version>[\w/_ -]+)", header)
        if rel_info is not None:
            rel_info = rel_info.group("version").decode()
        return rel_info

    def switch_collect_to_inference_phase(self):
        """method to switch from collect phase to inference phase.

        After calling this method, you can no longer call add_annotations
        on this object. From then on, only calls to add_inference are allowed."""
        # switch to inference phase (2). build index of annotations
        assert self.phase == 1
        self.flush()
        self.storage.create_index(self._anno_tab_name, ["EntryNr", "TermNr"])
        self.storage.__exit__(None, None, None)
        self.phase = 2
        # reopen in inference mode
        self.storage = TableStorageManager(self.inference_file, mode="w")
        self.storage.__enter__()
        self.storage.create_table(self.annotation_path, GeneOntologyTable, sort_key=["EntryNr", "TermNr", "Evidence"])

    def annotation_generated_date(self, date):
        self.go_tab.set_attr("annotations_generated", date)

    def add_annotations(self, enrs: Set[int], anno: GOA_Annotation):
        """parse go annotations and add them to the go buffer"""
        assert self.phase == 1, "cannot add annotations if phase != 1"
        if len(enrs) == 0:
            return
        try:
            term = self.go.term_by_id(anno.term_id)
        except ValueError:
            logger.warning(f"annotation {anno}: term {anno.term_id} not valid (obsolete?)")
            return

        ev = anno.evidence.encode("utf-8")
        ref = anno.db_ref.split("|")[0].encode("utf-8")
        data = [(enr, term.id, ev, ref) for enr in enrs]
        self.storage.add_rows(self._anno_tab_name, data)

    def add_inference(self, enr_term_tuples: Tuple[int, int], ref=b"OMA_Fun:001"):
        assert self.phase == 2, "cannot add inferences if phase != 2"
        data = []
        for enr, term_nr in enr_term_tuples:
            try:
                term = self.go.term_by_id(term_nr)
            except ValueError:
                logger.warning(f"annotation: term {term_nr} not valid (obsolete?)")
                continue
            data.append((enr, term.id, b"IEA", ref))
        self.storage.add_rows(self._anno_tab_name, data)


class XRefBasedMapper:
    def __init__(self, xref_db_path: os.PathLike):
        self.xref_db_path = xref_db_path

    def __enter__(self):
        self.xref_db = tables.open_file(self.xref_db_path, mode="r")
        self.xrefs: tables.Table = self.xref_db.get_node("/XRef")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.xref_db.close()

    @lru_cache(maxsize=65536)
    def find_id(self, id_):
        if m := re.match(r"(.*)\.\d{1,3}$", id_):
            id_ = m.group(1)
        query = id_.encode("utf-8")
        condvals = {"low": query, "high": query[:-1] + bytes([query[-1] + 1])}
        stmt = "(XRefId >= low) & (XRefId < high)"
        it = self.xrefs.where(stmt, condvals)
        xrefs = numpy.fromiter(
            (row.fetch_all_fields() for row in itertools.islice(it, 100)),
            dtype=self.xrefs.dtype,
        )
        if len(xrefs) == 0:
            return set([])
        xrefs.sort(order=["Verification", "EntryNr"])
        index = numpy.searchsorted(xrefs["Verification"], xrefs[0]["Verification"], side="right")
        return set(xrefs["EntryNr"][:index])


class Prof:
    __slots__ = ("t",)

    def __init__(self):
        self.t = {}

    def add(self, key, dt):
        self.t[key] = dt + self.t.get(key, 0.0)

    def __add__(self, other):
        """Support sum() over list of Prof; 0 + Prof = Prof"""
        if other == 0:
            return self
        new = self.__class__()
        for k, v in self.t.items():
            new.t[k] = v
        for k, v in other.t.items():
            new.t[k] = v + new.t.get(k, 0.0)
        return new

    def __str__(self):
        res = ["\nProfiling summary", "=" * 80]
        total = self.t.get("worker_total", 0.0)
        max_len = max(len(k) for k in self.t.keys())
        fmt_n = f"{{:{max_len}s}}: {{:12d}}"
        fmt_t = f"{{:{max_len}s}}: {{:10.2f}}s ({{:5.1f}}%)"
        for k, v in sorted(self.t.items()):
            if k.startswith("n_"):
                res.append(fmt_n.format(k, int(v)))
            else:
                pct = (v / total * 100) if total else 0
                res.append(fmt_t.format(k, v, pct))
        return "\n".join(res)


class Stats:
    def __init__(self):
        self.counts = {}

    def log(self, taxid, kind):
        if taxid not in self.counts:
            self.counts[taxid] = [0, 0, 0, 0]
        if kind < 0:
            self.counts[taxid][0] += 1
        elif kind == 0:
            self.counts[taxid][1] += 1
        elif kind == 1:
            self.counts[taxid][2] += 1
        if kind > 1:
            self.counts[taxid][3] += 1

    def __add__(self, other):
        """Return a new Stats object combining self + other"""
        new = self.__class__()
        # start with a copy of self.counts
        for taxid, counts in self.counts.items():
            new.counts[taxid] = counts.copy()
        # merge other
        for taxid, counts in other.counts.items():
            if taxid in new.counts:
                for i in range(4):
                    new.counts[taxid][i] += counts[i]
            else:
                new.counts[taxid] = counts.copy()
        return new

    def __radd__(self, other):
        """Support sum() over list of Stats; 0 + Stats = Stats"""
        if other == 0:
            return self
        return self.__add__(other)

    def summary(self):
        print("Mapping statistics")
        print("=" * 80)
        print(f"{'Taxon id':>15s}|{'NOT':>13s}|{'ID not found':>15s}|{'Found 1:1 match':>17s}|{'>1 ID matched':>15s}")
        print("-" * 80)
        for taxid in sorted(self.counts.keys(), key=lambda x: -sum(self.counts[x])):
            z = self.counts[taxid]
            print(f"{taxid:15d}|{z[0]:13d}|{z[1]:15d}|{z[2]:17d}|{z[3]:15d}")
        print("-" * 80)

    def __str__(self):
        return (
            f"Stats: 1:1: {sum(x[2] for x in self.counts.values())}; "
            f">1 ID: {sum(x[3] for x in self.counts.values())}; "
            f"NOT: {sum(x[0] for x in self.counts.values())}; "
            f"no-match: {sum(x[1] for x in self.counts.values())}; "
            f"[within {len(self.counts)} taxa]"
        )


def batch_gaf_annotations(gaf_files, batch_size=100000):
    batch = []
    for gaf in gaf_files:
        parser = AnnotationParser(gaf)
        for anno in parser:
            batch.append(anno)
            if len(batch) >= batch_size:
                yield batch
                batch = []
    if batch:
        yield batch


# module level globals for worker processes
_go_map_worker: Optional[XRefBasedMapper] = None
_relevant_taxid: Optional[Set[int]] = None
_fp_worker: Optional[FunctionPredictor] = None


def worker_logger_init(log_queue: Queue):
    # setup logging to use the queue
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers.clear()
    qh = QueueHandler(log_queue)
    root.addHandler(qh)


def _go_map_worker_init(xref_db_path, relevant_taxid: Set[int], log_queue: Queue):
    global _go_map_worker
    global _relevant_taxid
    worker_logger_init(log_queue)
    _go_map_worker = XRefBasedMapper(xref_db_path).__enter__()
    _relevant_taxid = relevant_taxid

    # Ensure __exit__ is called on exit
    Finalize(_go_map_worker, _go_map_worker.__exit__, args=(None, None, None), exitpriority=16)


def _go_map_worker_process(annotations):
    global _go_map_worker, _relevant_taxid

    stats = Stats()
    prof = Prof()
    results = []
    t0 = time.perf_counter()
    for annotation in annotations:
        t = time.perf_counter()
        # extract taxid (potentially contains also host taxon:1|taxon:1000)
        taxid = int(annotation.taxon.split("|")[0].split(":")[1])
        if taxid not in _relevant_taxid:
            continue
        if AnnotationFilter.is_negated(annotation):
            stats.log(taxid, -1)
            continue
        prof.add("preprocess_annotation", time.perf_counter() - t)

        t = time.perf_counter()
        enrs = _go_map_worker.find_id(annotation.db_obj_id)
        prof.add("find_xrefs", time.perf_counter() - t)

        stats.log(taxid, len(enrs))
        if not enrs:
            continue
        results.append((enrs, annotation))

    prof.add("worker_total", time.perf_counter() - t0)
    prof.add("n_annotations", len(annotations))
    logger.info(f"batch mapped: {stats}")
    return results, stats, prof


def _fp_worker_init(db_h5_path, go_h5_path, obo_path, clades: Optional[List[str]] = None, log_queue: Queue = None):
    global _fp_worker
    worker_logger_init(log_queue)

    gene_ontology = GeneOntology(OntologyParser(obo_path))
    gene_ontology.parse()

    _fp_worker = FunctionPredictor(
        db_h5_path=db_h5_path, go_h5_path=go_h5_path, ontology=gene_ontology, clades=clades
    ).__enter__()

    # Ensure __exit__ is called on exit
    Finalize(_fp_worker, _fp_worker.__exit__, args=(None, None, None), exitpriority=16)


def _fp_worker_process(oma_groups):
    global _fp_worker

    prof = Prof()
    results = []
    t0 = time.perf_counter()
    for group_nr in oma_groups:
        t = time.perf_counter()
        for enr, go_term in _fp_worker.annotate_group(group_nr):
            results.append((enr, go_term))
        prof.add("predict_group", time.perf_counter() - t)
    prof.add("worker_total", time.perf_counter() - t0)
    prof.add("n_groups", len(oma_groups))
    logger.info(f"batch predicted: {len(oma_groups)} groups, {len(results)} annotations")
    return results, prof


def stream_process_batches(
    pool,
    worker_func,
    batch_iter,
    max_inflight=2,
):
    """
    Submit batches lazily and yield completed results as soon as possible.

    max_inflight = number of outstanding futures allowed at once
    """
    inflight = collections.deque()

    for batch in batch_iter:
        # submit new batch
        inflight.append(pool.submit(worker_func, batch))

        # if too many inflight jobs, wait for one to finish
        if len(inflight) >= max_inflight:
            done, _ = wait(inflight, return_when=FIRST_COMPLETED)
            for fut in done:
                inflight.remove(fut)
                yield fut.result()

    # drain remaining futures
    while inflight:
        done, _ = wait(inflight, return_when=FIRST_COMPLETED)
        for fut in done:
            inflight.remove(fut)
            yield fut.result()


def setup_main_logging(log_queue):
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)-8s [pid=%(process)d name=%(processName)s] " "%(name)s: %(message)s"
    )
    handler.setFormatter(formatter)

    queue_listener = logging.handlers.QueueListener(log_queue, handler)
    queue_listener.start()
    return queue_listener


def import_go(
    obo: os.PathLike,
    gafs: List[os.PathLike],
    xref_db: os.PathLike,
    relevant_taxid: Set[int],
    og_db: os.PathLike,
    out: os.PathLike,
    clades: Optional[List[str]] = None,
    nr_procs: int = 1,
    batch_size: int = 100_000,
):
    stats = Stats()
    prof = Prof()

    log_queue = Queue()
    listener = setup_main_logging(log_queue)

    # attach to function for workers
    with GeneOntologyManager(out, "/Annotations/GeneOntology", "/Ontologies/GO") as go_man:
        go_man.add_ontology(obo_path=obo)
        with ProcessPoolExecutor(
            max_workers=nr_procs, initializer=_go_map_worker_init, initargs=(xref_db, relevant_taxid, log_queue)
        ) as pool:
            batch_iter = batch_gaf_annotations(gafs, batch_size=batch_size)

            for batch_results, batch_stats, c_prof in stream_process_batches(
                pool, _go_map_worker_process, batch_iter, max_inflight=nr_procs * 2
            ):
                stats += batch_stats
                prof += c_prof

                for enrs, anno in batch_results:
                    go_man.add_annotations(enrs, anno)

        stats.summary()
        go_man.switch_collect_to_inference_phase()

        with tables.open_file(og_db, "r") as db:
            tab: tables.Table = db.get_node("/Protein/Entries")
            idx = tab.colindexes["OmaGroup"][-1]
            nr_groups = int(tab[idx]["OmaGroup"])

        logger.info("predicting GO annotations based on OMA Groups")
        with ProcessPoolExecutor(
            max_workers=nr_procs,
            initializer=_fp_worker_init,
            initargs=(og_db, go_man.collect_file, obo, clades, log_queue),
        ) as pool:
            batch_iter = (range(i, min(i + 500, nr_groups + 1)) for i in range(1, nr_groups + 1, 500))
            for batch_results, c_prf in stream_process_batches(pool, _fp_worker_process, batch_iter):
                prof += c_prof
                go_man.add_inference(batch_results)

    logger.info("GO import completed")

    with DBBuilder(path=out, mode="append", logger=logger) as builder:
        builder.add_gene_ontology_term_cnts()
    print(prof)
