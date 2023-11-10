import itertools
import numpy
import collections
import tables
from .models import ProteinEntry


def chunked(iterator, dtype, n):
    """Batch data into numpy array of length n. The last batch may be shorter."""
    if n < 1:
        raise ValueError("n must be at least one")
    cnt_it, it = itertools.tee((row.fetch_all_fields() for row in iterator), 2)
    while True:
        counter = itertools.count()
        collections.deque(zip(itertools.islice(cnt_it, n), counter), maxlen=0)  # (consume at C speed)
        this_chunk_len = next(counter)
        if this_chunk_len == 0:
            break
        chunk = numpy.fromiter(it, dtype=dtype, count=this_chunk_len)
        yield chunk


def _count_annotations(h5, data):
    ex = numpy.array(
        [
            b"EXP",
            b"IDA",
            b"IPI",
            b"IMP",
            b"IGI",
            b"IEP",
            b"HTP",
            b"HDA",
            b"HMP",
            b"HGI",
            b"HEP",
        ],
        dtype="S3",
    )
    go_tab = h5.get_node("/Annotations/GeneOntology")
    it = go_tab.where('Evidence != b"IEA"')
    NG = len(data)
    for chunk in chunked(it, go_tab.dtype, 10000):
        genome_idx = data["EntryOff"].searchsorted(chunk["EntryNr"]) - 1
        ex_anno = numpy.where(numpy.isin(chunk["Evidence"], ex))[0]
        inc = numpy.bincount(genome_idx[ex_anno], minlength=NG)
        data["NrExpGO"] += inc

        cur_anno = numpy.where(numpy.isin(chunk["Evidence"], ex, invert=True))[0]
        inc = numpy.bincount(genome_idx[cur_anno], minlength=NG)
        data["NrCurGO"] += inc
    for pri, idx in enumerate(numpy.argsort(data, order=("NrExpGO", "NrCurGO"))[::-1]):
        data["Priority"][idx] = pri


class GenomeImportanceHelper:
    dtype = numpy.dtype([("EntryOff", "u4"), ("NrExpGO", "i4"), ("NrCurGO", "i4"), ("Priority", "i4")])
    tab_path = "/Summary/PerGenomeAnnotations"

    def __init__(self, h5: tables.File, compute_stats_if_needed=False):
        try:
            self.genome_data = h5.get_node(self.tab_path).read()
        except tables.NoSuchNodeError:
            self.genome_data = self._get_data(h5, compute_go=compute_stats_if_needed)

    @classmethod
    def create_and_store(cls, h5):
        inst = cls(h5, compute_stats_if_needed=True)
        try:
            tab = h5.get_node(cls.tab_path)
            h5.remove_node(tab)
        except tables.NoSuchNodeError:
            pass
        parent, name = cls.tab_path.rsplit("/", maxsplit=1)
        h5.create_table(
            parent,
            name,
            obj=inst.genome_data,
            createparents=True,
            expectedrows=len(inst.genome_data),
        )
        return inst

    def _get_data(self, h5, compute_go=True):
        offsets = h5.get_node("/Genome").read(field="EntryOff")
        data = numpy.zeros(len(offsets), dtype=self.dtype)
        data["EntryOff"] = offsets
        if compute_go:
            _count_annotations(h5, data)
        else:
            data["Priority"] = numpy.arange(len(data))
        return data

    def sort_numpy_with_entrynr(self, entries: numpy.ndarray):
        """This method sorts a numpy array with entry numbers (either column
        named EntryNr or a single integer array representing the entry numbers)
        according to the importance of the genome for experimentalists.

        :param entries: array containing entry numbers
        :returns: inplace sorted array
        """
        if not isinstance(entries, numpy.ndarray):
            raise ValueError("Invalid type: {}".format(type(entries)))

        if issubclass(entries.dtype.type, numpy.integer):
            col = entries
        elif "EntryNr" in entries.dtype.names:
            col = entries["EntryNr"]
        prios = self.genome_data[self.genome_data["EntryOff"].searchsorted(col) - 1]["Priority"]
        return entries[prios.argsort()]

    def sortkey(self, elements, entrynr_key=None):
        """get a sort key array that can be used to sort the elements
        according to the importance of the corresponding genome.

        :param elements: a list of elements, each corresponding to one protein
        :param entrynr_key: (optinal) a callable that extracts the entry_nr
            from the object. For ProteinModel instances as elements, the method
            automatically accesses the entry_nr using its property.
        """
        default_key = lambda x: x.entry_nr if isinstance(elements[0], ProteinEntry) else lambda x: x
        key = entrynr_key if entrynr_key is not None else default_key
        enrs = numpy.fromiter(map(key, elements), dtype=numpy.int32)
        prios = self.genome_data[self.genome_data["EntryOff"].searchsorted(enrs) - 1]["Priority"]
        return list(prios.argsort())
