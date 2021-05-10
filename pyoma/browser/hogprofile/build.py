import functools
import queue
import signal

import tables
import ete3
import multiprocessing as mp
import time
import tempfile
import pandas as pd
from datasketch import WeightedMinHashGenerator, MinHashLSHForest
from .pyhamutils import get_ham_treemap_from_row
from .hashutils import generate_treeweights, row2hash
from .tree_helper import (
    taxa_index,
    leaf_index,
    filter_tree,
    get_newick_tree_from_tax_db,
)
import logging
from logging.handlers import QueueHandler

logger = logging.getLogger(__name__)


class LSHBuilderBase(object):
    def __init__(self, db, numperm=256, taxfilter=None, taxmask=None, **kwargs):
        self.db = db
        tree_newick = get_newick_tree_from_tax_db(db.tax)
        full_tree = ete3.PhyloTree(tree_newick, format=1)
        self.tree = filter_tree(full_tree, root_node=taxmask, taxfilter=taxfilter)
        self.tree_newick = self.tree.write(format=1)
        self.taxa_index = taxa_index(self.tree)
        self.leaf_index = leaf_index(self.tree)
        self.tree_weights = generate_treeweights(self.tree, self.taxa_index)
        self.numperm = numperm
        self.wmg = WeightedMinHashGenerator(
            3 * len(self.taxa_index), sample_size=numperm, seed=1
        )


class BaseProfileBuilderProcess(mp.Process):
    def __init__(
        self,
        in_queue: mp.Queue,
        out_queue: mp.Queue,
        nr_workers_pre_step: mp.Value,
        log_queue: mp.Queue,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.in_queue = in_queue
        self.out_queue = out_queue
        self.log_queue = log_queue
        self.outstanding_dones = nr_workers_pre_step
        self.quit_req = False

    def setup(self):
        pass

    def finalize(self):
        pass

    def _handle_signal(self, signum, frame):
        print(
            "Stopping worker (pid: {}; name: {}): received signal {}".format(
                self.pid, self.name, signum
            )
        )
        self.quit_req = True

    def run(self):
        self.setup()
        # Set signals for worker process
        signal.signal(signal.SIGINT, self._handle_signal)
        signal.signal(signal.SIGUSR2, self._handle_signal)
        signal.signal(signal.SIGTERM, self._handle_signal)

        while not self.quit_req and self.outstanding_dones.value > 0:
            try:
                item = self.in_queue.get(timeout=1)
            except queue.Empty:
                continue

            if item is None:
                logger.info("received sentinel from previous pipeline step")
                with self.outstanding_dones.get_lock():
                    self.outstanding_dones.value -= 1
                logger.info("outstanding_dones: {}".format(self.outstanding_dones))

            else:
                result = self.handle_input(item)
                if result is not None:
                    self.out_queue.put(result)
        if not self.quit_req:
            self.out_queue.put(None)  # signal end of work for this worker
            logger.info("sent sentinel to next stage")
        else:
            logger.info("received termination signal")
        self.finalize()

    def handle_input(self, item):
        raise NotImplementedError("This method must be overwritten")


class SourceProcess(BaseProfileBuilderProcess):
    def __init__(self, out_queue, **kwargs):
        super().__init__(
            in_queue=None, out_queue=out_queue, nr_workers_pre_step=None, **kwargs
        )

    def _put_in_queue_with_status(self, input):
        while True:
            try:
                self.out_queue.put(input, timeout=5)
                break
            except queue.Full:
                qsize = "?"
                try:
                    qsize = self.out_queue.qsize()
                except NotImplementedError:
                    pass
                logger.info(
                    "Queue full. (approx size: {}) Cannot put item {} into output queue".format(
                        qsize, input
                    )
                )

    def run(self):
        self.setup()
        for input in self.generate_data():
            self._put_in_queue_with_status(input)
        self.out_queue.put(None)
        logger.info("sent sentinel to next stage")
        self.finalize()

    def generate_data(self):
        raise NotImplementedError("SourceProcess should implement this method")


class ProfileBuilder(BaseProfileBuilderProcess):
    def __init__(self, db_path, **kwargs):
        super().__init__(**kwargs)
        self.db_path = db_path
        self.builder = None
        self.ham_pipeline = None
        self.hash_pipeline = None

    def setup(self):
        from ..db import Database

        self.builder = LSHBuilderBase(db=Database(self.db_path))
        self.ham_pipeline = functools.partial(
            get_ham_treemap_from_row, tree=self.builder.tree_newick
        )
        self.hash_pipeline = functools.partial(
            row2hash,
            taxa_index=self.builder.taxa_index,
            species_index=self.builder.leaf_index,
            treeweights=self.builder.tree_weights,
            wmg=self.builder.wmg,
        )
        print("Worker {} init".format(self.name))

    def finalize(self):
        self.builder.db.close()

    def handle_input(self, df):
        print("handling df: {}".format(df))
        df["tree"] = df[["Fam", "ortho"]].apply(self.ham_pipeline, axis=1)
        df[["hash", "rows", "species"]] = df[["Fam", "tree"]].apply(
            self.hash_pipeline, axis=1
        )
        return df[["Fam", "hash", "species"]]


class Collector(BaseProfileBuilderProcess):
    def __init__(self, db_path, tmp_file, **kwargs):
        super().__init__(**kwargs)
        self.tmp_file = tmp_file
        self.db_path = db_path

    def _init_tmp_h5(self, root):
        print(root)
        hashes = self.h5.create_carray(
            root,
            "hashes",
            createparents=True,
            atom=tables.Int64Atom(),
            shape=(
                self.builder.db.get_nr_toplevel_hogs() + 1,
                self.builder.numperm * 2,
            ),
            filters=tables.Filters(complevel=3, complib="blosc", fletcher32=True),
        )
        species = self.h5.create_carray(
            root,
            "species_profile",
            createparents=True,
            atom=tables.Int8Atom(),
            shape=(
                self.builder.db.get_nr_toplevel_hogs() + 1,
                len(self.builder.leaf_index),
            ),
            filters=tables.Filters(complevel=3, complib="blosc", fletcher32=True),
        )
        lsh_forest = self.h5.create_vlarray(
            root,
            "min_hash_lsh_forest",
            createparents=True,
            atom=tables.ObjectAtom(),
            filters=tables.Filters(complevel=3, complib="blosc", fletcher32=True),
        )
        lsh_tree = self.h5.create_vlarray(
            root,
            "species_tree",
            createparents=True,
            atom=tables.ObjectAtom(),
            filters=tables.Filters(complevel=3, complib="blosc", fletcher32=True),
        )
        lsh_tree.append(self.builder.tree)
        self.h5.set_node_attr(root, "num_perm", self.builder.numperm)
        return hashes, species, lsh_forest

    def setup(self):
        from ..db import Database

        self.builder = LSHBuilderBase(db=Database(self.db_path))
        self.save_start = time.time()
        self.global_time = time.time()
        self.count = 0
        self.forest = MinHashLSHForest(num_perm=self.builder.numperm)
        self.h5 = tables.open_file(self.tmp_file, "w")
        root = "/HOGProfile/ALL"
        # if self.builder.tree.name != "0":
        #     root = "/HOGProfile/{}".format(self.builder.tree.name)
        print("storing results in {}:{}".format(self.tmp_file, root))

        self.hashes_matrix, self.species_matrix, self.forest_arr = self._init_tmp_h5(
            root
        )
        print("Collector process initialized")

    def handle_input(self, df: pd.DataFrame):
        print("handling hash df: {}".format(df))
        if not df.empty:
            hashes = df["hash"].to_dict()
            hashes = {fam: hashes[fam] for fam in hashes if hashes[fam]}
            for fam, hashvals in hashes.items():
                self.forest.add(str(fam), hashvals)
                self.hashes_matrix[fam, :] = hashvals.hashvalues.ravel()
                self.species_matrix[fam, :] = df["species"][fam]
                self.count += 1
            if time.time() - self.save_start > 200:
                print(
                    "Saving current results: {}".format(time.time() - self.global_time)
                )
                self.forest.index()
                print(
                    "Checking consistency of search index. Fam {} --> {}".format(
                        fam, self.forest.query(hashes[fam], k=10)
                    )
                )
                self.save_start = time.time()
                print("saved results")
        else:
            print("Empty dataframe in Collector: {}".format(df))

    def finalize(self):
        print("received all results. wrapping up...")
        print("computed minhashes: {}".format(self.count))
        self.forest.index()
        self.forest_arr.append(self.forest)
        self.h5.flush()
        self.h5.close()
        self.builder.db.close()
        print("stored all data in temporary file")


class HogGenerator(SourceProcess):
    def __init__(
        self, db_path, max_hogsize=None, min_hogsize=100, chunk_size=100, **kwargs
    ):
        super().__init__(**kwargs)
        self.db_path = db_path
        self.max_hogsize = max_hogsize
        self.min_hogsize = min_hogsize
        self.chunk_size = chunk_size

    def generate_data(self):
        from ..db import Database

        db = Database(self.db_path)
        nr_groups = db.get_nr_toplevel_hogs()
        chunk = {}
        logger.debug("nr of toplevel hogs: {}".format(nr_groups))
        for fam in range(1, nr_groups + 1):
            try:
                orthoxml = db.get_orthoxml(fam).decode()
            except ValueError as e:
                continue
            nr_species = orthoxml.count("<species name=")
            if (self.max_hogsize is None or nr_species < self.max_hogsize) and (
                self.min_hogsize is None or nr_species >= self.min_hogsize
            ):
                chunk[fam] = {"ortho": orthoxml}
            if len(chunk) > self.chunk_size:
                df = pd.DataFrame.from_dict(chunk, orient="index")
                df["Fam"] = df.index
                print("yielding chunk")
                yield df
                chunk = {}
        # send last chunk
        if len(chunk) > 0:
            df = pd.DataFrame.from_dict(chunk, orient="index")
            df["Fam"] = df.index
            yield df
        print("all hogs generated")
        db.close()


def rootlogger_configurer(queue):
    h = QueueHandler(queue)
    root = logging.getLogger()
    root.addHandler(h)
    root.setLevel(logging.DEBUG)


def logger_listener(log_queue):
    root = logging.getLogger()
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s"
    )
    console_handler.setFormatter(formatter)
    root.addHandler(console_handler)
    root.setLevel(logging.DEBUG)
    for record in iter(log_queue.get, None):
        logger = logging.getLogger(record.name)
        logger.handle(record)


class Stage(object):
    def __init__(self, process, nr_procs, **kwargs):
        self.process = process
        self.nr_procs = nr_procs
        self.kwargs = kwargs
        self.outstanding_dones = None
        self.in_queue = None
        self.out_queue = None


class Pipeline(object):
    def __init__(self):
        self.stages = []

    def add_stage(self, stage: Stage):
        self.stages.append(stage)

    def run(self):
        log_queue = mp.Queue()
        logger_process = mp.Process(target=logger_listener, args=(log_queue,))
        logger_process.start()
        rootlogger_configurer(log_queue)
        print("setting up pipeline")

        procs = []
        for i in range(len(self.stages) - 1, 0, -1):
            stage = self.stages[i]
            stage.outstanding_dones = mp.Value("d", self.stages[i - 1].nr_procs)
            stage.in_queue = mp.Queue(maxsize=10 * mp.cpu_count())
            stage.out_queue = (
                self.stages[i + 1].in_queue if i < len(self.stages) - 1 else mp.Queue()
            )
            for k in range(stage.nr_procs):
                p = stage.process(
                    in_queue=stage.in_queue,
                    out_queue=stage.out_queue,
                    nr_workers_pre_step=stage.outstanding_dones,
                    log_queue=log_queue,
                    daemon=True,
                    **stage.kwargs
                )
                procs.append(p)
                p.start()
        # add source process
        stage = self.stages[0]
        stage.out_queue = self.stages[1].in_queue
        for k in range(stage.nr_procs):
            p = stage.process(
                out_queue=stage.out_queue, log_queue=log_queue, **stage.kwargs
            )
            procs.append(p)
            p.start()
        print("all processes started")
        try:
            for p in procs:
                p.join()
        except KeyboardInterrupt:
            print("keyboard interupt in main loop")
            time.sleep(20)
        log_queue.put(None)
        logger_process.join()

        print("successfully joined all processes")


def compute_profiles(db_path, min_hogsize=100, max_hogsize=None, nr_procs=None):
    pipeline = Pipeline()
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as h5_tmp:
        tmp_file = h5_tmp.name
    if nr_procs is None:
        nr_procs = mp.cpu_count()

    pipeline.add_stage(
        Stage(
            HogGenerator,
            nr_procs=1,
            db_path=db_path,
            min_hogsize=min_hogsize,
            max_hogsize=max_hogsize,
        )
    )
    pipeline.add_stage(Stage(ProfileBuilder, nr_procs=nr_procs, db_path=db_path))
    pipeline.add_stage(Stage(Collector, nr_procs=1, db_path=db_path, tmp_file=tmp_file))
    print("generated pipeline. about to starting it")
    pipeline.run()
    print("finished computing profiles")

    with tables.open_file(db_path, "a") as db, tables.open_file(tmp_file, "r") as tmp:
        tmp.root._f_copy_children(db.root, recursive=True, overwrite=True)
    print("Finished writing everything")
