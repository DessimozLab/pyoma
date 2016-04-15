from __future__ import division
from future.builtins import str
from future.builtins import chr
from future.builtins import range
from future.builtins import object
from future.builtins import super
from future.standard_library import hooks
import csv
import resource
import tables
import numpy
import numpy.lib.recfunctions
import os
import subprocess
import errno
import json
import time
import familyanalyzer
import re
import multiprocessing as mp
import lxml.html
import collections
import gzip
import hashlib
import itertools

from .. import common
from . import locus_parser
from . import tablefmt
with hooks():
    import urllib.request


class DarwinException(Exception):
    pass


def callDarwinExport(func, drwfile=None):
    """Function starts a darwin session, loads convert.drw file
    and calls the darwin function passed as argument. The output
    is expected to be written by darwin in json format into the
    file specified by 'outfn'.
    This function returns the parsed json datastructure"""

    tmpfile = "/tmp/darwinExporter_{:d}.dat".format(os.getpid())
    if drwfile is None:
        drwfile = os.path.abspath(os.path.splitext(__file__)[0] + ".drw")
    try:
        with open(os.devnull, 'w') as DEVNULL:
            stacksize = resource.getrlimit(resource.RLIMIT_STACK)
            common.package_logger.info('current stacklimit: {}'.format(stacksize))
            common.package_logger.info('setting stacklimit: {}'.format((max(stacksize)-1, stacksize[1])))
            resource.setrlimit(resource.RLIMIT_STACK, (min(stacksize), stacksize[1]))
            p = subprocess.Popen(['darwin', '-q', '-E', '-B'], stdin=subprocess.PIPE,
                                 stderr=subprocess.PIPE, stdout=DEVNULL)
            drw_cmd = "outfn := '{}': ReadProgram('{}'): {}; done;".format(
                tmpfile, drwfile, func).encode('utf-8')
            common.package_logger.debug('calling darwin function: {}'.format(func))
            p.communicate(input=drw_cmd)
            if p.returncode > 0:
                raise DarwinException(p.stderr.read())

        trans_tab = "".join(str(chr(x)) for x in range(128)) + " " * 128
        with open(tmpfile, 'r') as jsonData:
            rawdata = jsonData.read()
            data = json.loads(rawdata.translate(trans_tab))
    finally:
        silentremove(tmpfile)
    return data


def uniq(seq):
    """return uniq elements of a list, preserving order

    :param seq: an iterable to be analyzed
    """
    seen = set()
    return [x for x in seq if not (x in seen or seen.add(x))]


def silentremove(filename):
    """Function to remove a given file. No exception is raised if the
    file does not exist. Other errors are passed to the user.
    :param filename: the path of the file to be removed"""
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured


def load_tsv_to_numpy(args):
    fn, off1, off2, swap = args
    relEnum = tablefmt.PairwiseRelationTable.columns['RelType'].enum._names
    relEnum['n:1'] = relEnum['m:1']
    relEnum['1:m'] = relEnum['1:n']
    relEnum['n:m'] = relEnum['m:n']
    read_dir = -1 if swap else 1
    dtype = [('EntryNr1', 'i4'), ('EntryNr2', 'i4'), ('Score', 'f4'), ('RelType', 'i1'),
             ('AlignmentOverlap', 'f2'), ('Distance', 'f4')]
    for curNr, curFn in enumerate([fn, fn.replace('.ext.', '.')]):
        try:
            with gzip.GzipFile(curFn) as fh:
                augData = numpy.genfromtxt(fh, dtype=dtype,
                                           names=[_[0] for _ in dtype],
                                           delimiter='\t',
                                           usecols=(0, 1, 2, 3, 4, 5),
                                           converters={'EntryNr1': lambda nr: int(nr) + off1,
                                                       'EntryNr2': lambda nr: int(nr) + off2,
                                                       'RelType': lambda rel: (relEnum[rel[::read_dir].decode()]
                                                                               if len(rel) <= 3
                                                                               else relEnum[rel.decode()]),
                                                       'Score': lambda score: float(score)/100})
                break
        except OSError as e:
            if curNr < 1:
                common.package_logger.info('tried to load {}'.format(curFn))
                pass
            else:
                raise e

    if swap:
        reversed_cols = tuple(augData.dtype.names[z] for z in (1, 0, 2, 3, 4, 5))
        augData.dtype.names = reversed_cols
    full_table = numpy.empty(augData.size, dtype=tables.dtype_from_descr(tablefmt.PairwiseRelationTable))
    full_table[0:] = augData
    for col_not_in_tsv in set(full_table.dtype.names) - set(augData.dtype.names):
        full_table[col_not_in_tsv] = -1
    return full_table


def read_vps_from_tsv(gs, ref_genome):
    ref_genome_idx = gs.get_where_list('(UniProtSpeciesCode=={!r})'.
                                       format(ref_genome))[0]
    job_args = []
    for g in range(len(gs)):
        if g == ref_genome_idx:
            continue
        g1, g2 = sorted((g, ref_genome_idx,))
        off1, off2 = gs.read_coordinates(numpy.array((g1, g2)), 'EntryOff')
        fn = os.path.join(os.environ['DARWIN_OMADATA_PATH'], 'Phase4',
                          gs.cols.UniProtSpeciesCode[g1].decode(),
                          gs.cols.UniProtSpeciesCode[g2].decode() + ".orth.txt.gz")
        tup = (fn, off1, off2, g1 != ref_genome_idx)
        common.package_logger.info('adding job: {}'.format(tup))
        job_args.append(tup)

    pool = mp.Pool()
    all_pairs = pool.map(load_tsv_to_numpy, job_args)
    pool.close()
    return numpy.lib.recfunctions.stack_arrays(all_pairs, usemask=False)


class DataImportError(Exception):
    pass


class DarwinExporter(object):
    DB_SCHEMA_VERSION = '2.0'
    DRW_CONVERT_FILE = os.path.abspath(os.path.splitext(__file__)[0] + '.drw')

    def __init__(self, path, logger=None):
        self.logger = logger if logger is not None else common.package_logger
        fn = os.path.normpath(os.path.join(
            os.environ['DARWIN_BROWSERDATA_PATH'],
            path))
        mode = 'append' if os.path.exists(fn) else 'write'
        self._compr = tables.Filters(complevel=6, complib='zlib', fletcher32=True)
        self.h5 = tables.open_file(fn, mode=mode[0], filters=self._compr)
        self.logger.info("opened {} in {} mode, options {}".format(
            fn, mode, str(self._compr)))
        if mode == 'write':
            self.h5.root._f_setattr('convertion_start', time.strftime("%c"))

    def call_darwin_export(self, func):
        return callDarwinExport(func, self.DRW_CONVERT_FILE)

    def _get_or_create_node(self, path, desc=None):
        try:
            grp = self.h5.get_node(path)
        except tables.NoSuchNodeError:
            base, name = os.path.split(path)
            grp = self.h5.create_group(base, name, title=desc, createparents=True)
        return grp

    def get_version(self):
        """return version of the dataset.

        Default implementation searches for 'mname' in Matrix or matrix_stats.drw files.
        """
        for fname in ('Matrix', 'matrix_stats.drw'):
            with open(os.path.join(os.environ['DARWIN_BROWSERDATA_PATH'], fname), 'r') as fh:
                for i, line in enumerate(fh):
                    if line.startswith('mname :='):
                        match = re.match(r'mname := \'(?P<version>[^\']*)\'', line)
                        return match.group('version')
                    if i > 1000:
                        break
        raise DataImportError('No version information found')

    def add_version(self):
        version = self.get_version()
        self.h5.set_node_attr('/', 'oma_version', version)
        self.h5.set_node_attr('/', 'pytables', tables.get_pytables_version())
        self.h5.set_node_attr('/', 'hdf5_version', tables.get_hdf5_version())
        self.h5.set_node_attr('/', 'db_schema_version', self.DB_SCHEMA_VERSION)

    def add_species_data(self):
        cache_file = os.path.join(
            os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
            'pyoma', 'gs.json')
        if os.path.exists(cache_file):
            with open(cache_file, 'r') as fd:
                data = json.load(fd)
        else:
            data = self.call_darwin_export('GetGenomeData();')
        gstab = self.h5.create_table('/', 'Genome', tablefmt.GenomeTable,
                                     expectedrows=len(data['GS']))
        self._write_to_table(gstab, data['GS'])
        gstab.cols.NCBITaxonId.create_csindex(filters=self._compr)
        gstab.cols.UniProtSpeciesCode.create_csindex(filters=self._compr)
        gstab.cols.EntryOff.create_csindex(filters=self._compr)

        taxtab = self.h5.create_table('/', 'Taxonomy', tablefmt.TaxonomyTable,
                                      expectedrows=len(data['Tax']))
        self._write_to_table(taxtab, data['Tax'])
        taxtab.cols.NCBITaxonId.create_csindex(filters=self._compr)

    def _convert_to_numpyarray(self, data, tab):
        """convert a list of list dataset into a numpy rec array that
        corresponds to the table definition of `tab`.

        :param data: the data to be converted.
        :param tab: a pytables table node."""

        enum_cols = {i: tab.get_enum(col) for (i, col) in enumerate(tab.colnames)
                     if tab.coltypes[col] == 'enum'}
        dflts = [tab.coldflts[col] for col in tab.colnames]

        def map_data(col, data):
            try:
                val = data[col]
                return enum_cols[col][val]
            except IndexError:
                return dflts[col]
            except KeyError:
                return val

        arr = numpy.empty(len(data), dtype=tab.dtype)
        for i, row in enumerate(data):
            as_tup = tuple(map_data(c, row) for c in range(len(dflts)))
            arr[i] = as_tup
        return arr

    def add_orthologs(self):
        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            rel_node_for_genome = self._get_or_create_node('/PairwiseRelation/{}'.format(genome))
            if 'VPairs' not in rel_node_for_genome:
                cache_file = os.path.join(
                    os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
                    'pyoma', 'vps', '{}.json'.format(genome))
                if os.path.exists(cache_file):
                    with open(cache_file, 'r') as fd:
                        data = json.load(fd)
                elif ((not os.getenv('DARWIN_OMADATA_PATH') is None) and
                          os.path.exists(os.path.join(
                                  os.environ['DARWIN_OMADATA_PATH'], 'Phase4'))):
                    # try to read from Phase4 in parallel.
                    data = read_vps_from_tsv(self.h5.root.Genome,
                                             genome.encode('utf-8'))
                else:
                    # fallback to read from VPsDB
                    data = self.call_darwin_export('GetVPsForGenome({})'.format(genome))

                vp_tab = self.h5.create_table(rel_node_for_genome, 'VPairs', tablefmt.PairwiseRelationTable,
                                              expectedrows=len(data))
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, vp_tab)
                self._write_to_table(vp_tab, data)
                vp_tab.cols.EntryNr1.create_csindex()

    def add_same_species_relations(self):
        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            rel_node_for_genome = self._get_or_create_node('/PairwiseRelation/{}'.format(genome))
            if 'within' not in rel_node_for_genome:
                cache_file = os.path.join(
                    os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
                    'pyoma', 'cps', '{}.json'.format(genome))
                if os.path.exists(cache_file):
                    with open(cache_file, 'r') as fd:
                        data = json.load(fd)
                else:
                    # fallback to read from VPsDB
                    data = self.call_darwin_export('GetSameSpeciesRelations({})'.format(genome))

                ss_tab = self.h5.create_table(rel_node_for_genome, 'within', tablefmt.PairwiseRelationTable,
                                              expectedrows=len(data))
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, ss_tab)
                self._write_to_table(ss_tab, data)
                ss_tab.cols.EntryNr1.create_csindex()

    def _add_sequence(self, sequence, row, sequence_array, off, typ="Seq"):
        # add ' ' after each sequence (Ascii is smaller than
        # any AA, allows to build PAT array with split between
        # sequences.
        seqLen = len(sequence) + 1
        row[typ + 'BufferOffset'] = off
        row[typ + 'BufferLength'] = seqLen
        seqNumpyObj = numpy.ndarray((seqLen,),
                                    buffer=(sequence + " ").encode('utf-8'),
                                    dtype=tables.StringAtom(1))
        sequence_array.append(seqNumpyObj)
        if typ == "Seq":
            row['MD5ProteinHash'] = hashlib.md5(sequence.encode('utf-8')).hexdigest()
        return seqLen

    def add_proteins(self):
        gsNode = self.h5.get_node('/Genome')
        nrProt = sum(gsNode.cols.TotEntries)
        nrAA = sum(gsNode.cols.TotAA)
        try:
            protGrp = self.h5.get_node('/Protein')
        except tables.NoSuchNodeError:
            protGrp = self.h5.create_group('/', 'Protein')
        protTab = self.h5.create_table(protGrp, 'Entries', tablefmt.ProteinTable,
                                       expectedrows=nrProt)
        seqArr = self.h5.create_earray(protGrp, 'SequenceBuffer',
                                       tables.StringAtom(1), (0,), 'concatenated protein sequences',
                                       expectedrows=nrAA + nrProt)
        cdnaArr = self.h5.create_earray(protGrp, 'CDNABuffer',
                                        tables.StringAtom(1), (0,), 'concatenated cDNA sequences',
                                        expectedrows=3 * nrAA + nrProt)
        seqOff = cdnaOff = 0
        for gs in gsNode.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            cache_file = os.path.join(
                os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
                'pyoma', 'prots', '{}.json'.format(genome))
            if os.path.exists(cache_file):
                with open(cache_file, 'r') as fd:
                    data = json.load(fd)
            else:
                data = self.call_darwin_export('GetProteinsForGenome({})'.format(genome))

            if len(data['seqs']) != gs['TotEntries']:
                raise DataImportError('number of entries ({:d}) does '
                                      'not match number of seqs ({:d}) for {}'.format(
                                        len(data['seqs']), gs['TotEntries'], genome))

            locTab = self.h5.create_table('/Protein/Locus',
                                          genome, tablefmt.LocusTable, createparents=True,
                                          expectedrows=gs['TotEntries'] * 4)

            for nr in range(gs['TotEntries']):
                eNr = data['off'] + nr + 1
                protTab.row['EntryNr'] = eNr
                protTab.row['OmaGroup'] = data['ogs'][nr]

                seqOff += self._add_sequence(data['seqs'][nr], protTab.row, seqArr, seqOff)
                cdnaOff += self._add_sequence(data['cdna'][nr], protTab.row, cdnaArr, cdnaOff, 'CDNA')

                protTab.row['Chromosome'] = data['chrs'][nr]
                protTab.row['AltSpliceVariant'] = data['alts'][nr]
                protTab.row['OmaHOG'] = b" "  # will be assigned later
                protTab.row['CanonicalId'] = b" "  # will be assigned later

                locus_str = data['locs'][nr]
                try:
                    locus_tab = locus_parser.parse(eNr, locus_str)
                    locTab.append(locus_tab)

                    protTab.row['LocusStart'] = locus_tab.Start.min()
                    protTab.row['LocusEnd'] = locus_tab.End.max()
                    protTab.row['LocusStrand'] = locus_tab[0].Strand
                except ValueError as e:
                    self.logger.warning(e)
                protTab.row['SubGenome'] = data['subgenome'][nr].encode('ascii')
                protTab.row.append()
            protTab.flush()
            seqArr.flush()
            for n in (protTab, seqArr, locTab):
                self.logger.info('worte %s: compression ratio %3f%%' %
                                 (n._v_pathname, 100 * n.size_on_disk / n.size_in_memory))
        protTab.cols.EntryNr.create_csindex(filters=self._compr)
        protTab.cols.MD5ProteinHash.create_csindex(filters=self._compr)

    def _write_to_table(self, tab, data):
        tab.append(data)
        self.logger.info('wrote %s : compression ratio %.3f%%' %
                         (tab._v_pathname, 100 * tab.size_on_disk / tab.size_in_memory))

    def add_hogs(self):
        hog_path = os.path.normpath(os.path.join(
            os.environ['DARWIN_BROWSERDATA_PATH'],
            '..', 'downloads', 'HOGs'))
        entryTab = self.h5.get_node('/Protein/Entries')
        tree_filename = os.path.join(
            os.environ['DARWIN_BROWSERDATA_PATH'],
            'speciestree.nwk')
        hog_converter = HogConverter(entryTab)
        hog_converter.attach_newick_taxonomy(tree_filename)
        hogTab = self.h5.create_table('/', 'HogLevel', tablefmt.HOGsTable,
                                      'nesting structure for each HOG', expectedrows=1e8)
        self.orthoxml_buffer = self.h5.create_earray('/OrthoXML', 'Buffer',
                                                     tables.StringAtom(1), (0,), 'concatenated orthoxml files',
                                                     expectedrows=1e9, createparents=True)
        self.orthoxml_index = self.h5.create_table('/OrthoXML', 'Index', tablefmt.OrthoXmlHogTable,
                                                   'Range index per HOG into OrthoXML Buffer', expectedrows=5e6)
        for root, dirs, filenames in os.walk(hog_path):
            for fn in filenames:
                try:
                    levels = hog_converter.convert_file(os.path.join(root, fn))
                    hogTab.append(levels)
                    fam_nrs = set([z[0] for z in levels])
                    self.add_orthoxml(os.path.join(root, fn), fam_nrs)
                except Exception as e:
                    self.logger.error('an error occured while processing ' + fn + ':')
                    self.logger.exception(e)

        hog_converter.write_hogs()

    def add_orthoxml(self, orthoxml_path, fam_nrs):
        """append orthoxml file content to orthoxml_buffer array and add index for the HOG family"""
        if len(fam_nrs) > 1:
            self.logger.warning('expected only one family per HOG, but found {}: {}\n'.format(len(fam_nrs), fam_nrs))
        with open(orthoxml_path, 'r') as fh:
            orthoxml = fh.read().encode('utf-8')
            offset = len(self.orthoxml_buffer)
            length = len(orthoxml)

            self.orthoxml_buffer.append(numpy.ndarray((length,),
                buffer=orthoxml, dtype=tables.StringAtom(1)))
            for fam in fam_nrs:
                row = self.orthoxml_index.row
                row['Fam'] = fam
                row['HogBufferOffset'] = offset
                row['HogBufferLength'] = length
                offset += length
                row.append()

    def add_xrefs(self):
        self.logger.info('start parsing ServerIndexed to extract XRefs and GO annotations')
        db_parser = IndexedServerParser()
        xref_tab = self.h5.create_table('/', 'XRef', tablefmt.XRefTable,
                                        'Cross-references of proteins to external ids / descriptions',
                                        expectedrows=1e8)
        go_tab = self.h5.create_table('/', 'GeneOntology', tablefmt.GeneOntologyTable,
                                      'Gene Ontology annotations', expectedrows=1e8)
        with DescriptionManager(self.h5, '/Protein/Entries', '/Protein/DescriptionBuffer') as de_man:
            xref_importer = XRefImporter(db_parser, xref_tab, go_tab, de_man)
            db_parser.parse_entrytags()
            xref_importer.flush_buffers()

    def close(self):
        self.h5.root._f_setattr('conversion_end', time.strftime("%c"))
        self.h5.close()
        self.logger.info('closed {}'.format(self.h5.filename))

    def create_indexes(self):
        self.logger.info('createing indexes for HOGs')
        hogTab = self.h5.get_node('/HogLevel')
        hogTab.cols.Fam.create_index()
        hogTab.cols.ID.create_index()
        orthoxmlTab = self.h5.get_node('/OrthoXML/Index')
        orthoxmlTab.cols.Fam.create_csindex()
        entryTab = self.h5.get_node('/Protein/Entries')
        entryTab.cols.OmaHOG.create_csindex()

        self.logger.info('creating index for xrefs (EntryNr and XRefId)')
        xrefTab = self.h5.get_node('/XRef')
        xrefTab.cols.EntryNr.create_csindex()
        xrefTab.cols.XRefId.create_csindex()

        self.logger.info('creating index for go (EntryNr and TermNr)')
        goTab = self.h5.get_node('/GeneOntology')
        goTab.cols.EntryNr.create_csindex()
        goTab.cols.TermNr.create_index()

        self.logger.info('creating index for domains (EntryNr)')
        domtab = self.h5.get_node('/Annotations/Domains')
        domtab.cols.EntryNr.create_csindex()

    def _iter_canonical_xref(self):
        """extract one canonical xref id for each protein.

        We take the first valid xref per gene with the ordering of xrefsources
        as given in the xrefsource_order."""
        xrefsource_order = ('UniProtKB/SwissProt', 'UniProtKB/TrEMBL',
                            'Ensembl Gene', 'Ensembl Protein', 'FlyBase',)

        xrefs = self.h5.get_node('/XRef')
        source_enum = xrefs.get_enum('XRefSource')
        canonical_sources = [source_enum[z] for z in xrefsource_order]
        current_protein = None
        past_proteins = set([])
        for xref in xrefs:
            if xref['EntryNr'] != current_protein:
                if current_protein:
                    past_proteins.add(current_protein)
                    yield (current_protein, current_xrefs.pop()[1])
                current_protein = xref['EntryNr']
                current_xrefs = [(1000, b'')]  # init with a sentinel
                if current_protein in past_proteins:
                    raise DataImportError('Data in /XRef is not sorted w.r.t. EntryNr')
            if xref['XRefSource'] in canonical_sources and xref['XRefSource'] < current_xrefs[-1][0]:
                current_xrefs.append((xref['XRefSource'], xref['XRefId']))
        if current_protein:
            yield (current_protein, current_xrefs.pop()[1])

    def add_canonical_id(self):
        """add one canonical xref id to the /Protein/Entries table."""
        self.logger.info('adding canonical ids for each protein...')
        prot_tab = self.h5.get_node('/Protein/Entries')
        canonical_ids = numpy.chararray(shape=(len(prot_tab),), itemsize=prot_tab.cols.CanonicalId.dtype.itemsize)
        for eNr, canonical_id in self._iter_canonical_xref():
            row_nr = eNr - 1
            row = prot_tab[row_nr]
            if row['EntryNr'] != eNr:
                self.logger.warn('Entries table not properly sorted: {}, expected {}'.format(row['EntryNr'], eNr))
                raise DataImportError('Entries table not properly sorted')
            canonical_ids[row_nr] = canonical_id
        prot_tab.modify_column(0, len(prot_tab), 1, column=canonical_ids, colname='CanonicalId')
        prot_tab.flush()

    def add_domain_info(self, domains):
        self.logger.info('adding domain information...')
        domtab = self.h5.create_table('/Annotations', 'Domains', tablefmt.DomainTable, createparents=True, expectedrows=1e7)
        entrytab = self.h5.get_node('/Protein/Entries')
        md5_to_enr = collections.defaultdict(list)
        for e in entrytab:
            md5_to_enr[e['MD5ProteinHash']].append(e['EntryNr'])

        buffer = []
        for i, domain in enumerate(domains):
            for entry_nr in md5_to_enr[domain.md5.encode('utf-8')]:
                buffer.append((entry_nr, domain.id, domain.coords))
                if len(buffer) > 5000:
                    domtab.append(buffer)
                    buffer = []
            if i % 50000 == 0:
                self.logger.info('processed {:d} domain annotations so far'.format(i))
        if len(buffer) > 0:
            domtab.append(buffer)
        domtab.flush()

    def add_domainname_info(self, domainname_infos):
        self.logger.info('adding domain name information...')
        dom_name_tab = self.h5.create_table('/Annotations', 'DomainDescription', tablefmt.DomainDescriptionTable,
                                            createparents=True, expectedrows=2e5)
        buffer = []
        for i, dom_info in enumerate(domainname_infos):
            buffer.append(dom_info)
            if len(buffer) > 5000:
                self._write_to_table(dom_name_tab, buffer)
                buffer = []
            if i % 50000 == 0:
                self.logger.info('processed {:d} domain name descriptions so far'.format(i))
        if len(buffer) > 0:
            self._write_to_table(dom_name_tab, buffer)
        dom_name_tab.flush()


def download_url_if_not_present(url):
    fname = os.path.join(os.getenv('DARWIN_NETWORK_SCRATCH_PATH', '/tmp'), "Browser", "xref",
                         url.split('/')[-1])
    if not os.path.exists(fname):
        try:
            urllib.request.urlretrieve(url, fname)
        except urllib.request.URLError:
            common.package_logger.warn('cannot access domain url')
    return fname


def iter_domains(url):
    DomainTuple = collections.namedtuple('DomainTuple', ('md5', 'id', 'coords'))

    fname = download_url_if_not_present(url)
    with gzip.open(fname, 'rt') as uncompressed:
        csv_reader = csv.reader(uncompressed)
        for lineNr, row in enumerate(csv_reader):
            try:
                dom = DomainTuple(*row[0:3])
                yield dom
            except Exception:
                common.package_logger.exception('cannot create tuple from line {}'.format(lineNr))


def only_pfam_or_cath_domains(iterable):
    cath_re = re.compile(r'[1-4]\.')
    for dom in iterable:
        if dom.id.startswith('PF') or cath_re.match(dom.id) is not None:
            yield dom


class DescriptionManager(object):
    def __init__(self, db, entry_path, buffer_path):
        self.db = db
        self.entry_path = entry_path
        self.buffer_path = buffer_path

    def __enter__(self):
        self.entry_tab = self.db.get_node(self.entry_path)
        if not numpy.all(numpy.equal(self.entry_tab.col('EntryNr'),
                                     numpy.arange(1,len(self.entry_tab)+1))):
            raise RuntimeError('entry table is not sorted')

        root, name = os.path.split(self.buffer_path)
        self.desc_buf = self.db.create_earray(root, name,
            tables.StringAtom(1), (0,), 'concatenated protein descriptions',
            expectedrows=len(self.entry_tab)*100)
        self.cur_eNr = None
        self.cur_desc = []
        bufindex_dtype = numpy.dtype([(col, self.entry_tab.coldtypes[col])
            for col in ('DescriptionOffset', 'DescriptionLength')])
        # columns to be stored in entry table with buffer index data
        self.buf_index = numpy.zeros(len(self.entry_tab), dtype=bufindex_dtype)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.cur_eNr:
            self._store_description()
        self.desc_buf.flush()
        self.entry_tab.modify_columns(columns=self.buf_index,
                                      names=self.buf_index.dtype.names)
        self.entry_tab.flush()

    def add_description(self, eNr, desc):
        """stages a description for addition. Note that the descriptions 
        must be ordered according to the entryNr, i.e. all descriptions 
        related to eNr X must be staged before changeing to another eNr."""
        if self.cur_eNr and self.cur_eNr != eNr:
            self._store_description()
            self.cur_desc = []
        self.cur_eNr = eNr
        self.cur_desc.append(desc)

    def _store_description(self):
        buf = "; ".join(self.cur_desc).encode('utf-8')
        buf = buf[0:2**16-1]  # limit to max value of buffer length field
        len_buf = len(buf)
        idx = self.cur_eNr - 1
        self.buf_index[idx]['DescriptionOffset'] = len(self.desc_buf)
        self.buf_index[idx]['DescriptionLength'] = len_buf
        self.desc_buf.append(numpy.ndarray((len_buf,), buffer=buf, dtype=tables.StringAtom(1)))


class GroupAnnotatorInclGeneRefs(familyanalyzer.GroupAnnotator):
    def _annotateGroupR(self, node, og, idx=0):
        if familyanalyzer.OrthoXMLQuery.is_geneRef_node(node):
            node.set('og', og)
        else:
            super()._annotateGroupR(node, og, idx)


class HogConverter(object):
    def __init__(self, entry_tab):
        self.fam_re = re.compile(r'HOG:(?P<fam_nr>\d+)')
        self.hogs = numpy.zeros(shape=(len(entry_tab) + 1,), dtype=entry_tab.cols.OmaHOG.dtype)
        self.entry_tab = entry_tab

    def attach_newick_taxonomy(self, tree):
        self.taxonomy = familyanalyzer.NewickTaxonomy(tree)

    def convert_file(self, fn):
        p = familyanalyzer.OrthoXMLParser(fn)
        if self.taxonomy:
            p.augmentTaxonomyInfo(self.taxonomy)
        GroupAnnotatorInclGeneRefs(p).annotateDoc()

        levs = []
        for fam in p.getToplevelGroups():
            m = self.fam_re.match(fam.get('og'))
            fam_nr = int(m.group('fam_nr'))
            levs.extend([(fam_nr, n.getparent().get('og'), n.get('value'),)
                         for n in p._findSubNodes('property')
                         if n.get('name') == "TaxRange"])

        geneNodes = p.root.findall('.//{{{ns0}}}geneRef'.
                                   format(**familyanalyzer.OrthoXMLParser.ns))
        for x in geneNodes:
            self.hogs[int(x.get('id'))] = x.get('og')

        return levs

    def write_hogs(self):
        self.entry_tab.modify_column(0, len(self.entry_tab), 1, self.hogs[1:], 'OmaHOG')
        self.entry_tab.flush()


class XRefImporter(object):
    def __init__(self, db_parser, xref_tab, go_tab, desc_manager):
        self.xrefs = []
        self.go = []
        self.xref_tab = xref_tab
        self.go_tab = go_tab
        self.desc_manager = desc_manager

        xrefEnum = tablefmt.XRefTable.columns.get('XRefSource').enum
        for tag in ['GI', 'EntrezGene', 'WikiGene', 'IPI']:
            db_parser.add_tag_handler(
                tag,
                lambda key, enr, typ=xrefEnum[tag]: self.multi_key_handler(key, enr, typ))
        db_parser.add_tag_handler('Refseq_ID',
                                  lambda key, enr: self.multi_key_handler(key, enr, xrefEnum['RefSeq']))
        db_parser.add_tag_handler('SwissProt',
                                  lambda key, enr: self.key_value_handler(key, enr, xrefEnum['UniProtKB/SwissProt']))
        db_parser.add_tag_handler('DE',
                                  lambda key, enr: self.description_handler(key, enr))
        db_parser.add_tag_handler('GO',self.go_handler)
        db_parser.add_tag_handler('ID', self.assign_source_handler)
        db_parser.add_tag_handler('AC', self.assign_source_handler)

        db_parser.add_tag_handler('ID',
                                  lambda key, enr: self.multi_key_handler(key, enr, xrefEnum['SourceID']))
        db_parser.add_tag_handler('AC',
                                  lambda key, enr: self.multi_key_handler(key, enr, xrefEnum['SourceAC']))
        for tag in ['PMP', 'EMBL']:
            db_parser.add_tag_handler(
                tag,
                lambda key, enr, typ=xrefEnum[tag]: self.multi_key_handler(key, enr, typ))
        for tag in ['SwissProt_AC', 'UniProt']:  # UniProt/TrEMBL tag is cut to UniProt!
            db_parser.add_tag_handler(tag,
                                      lambda key, enr, typ=xrefEnum['UniProtKB/TrEMBL']:
                                      self.remove_uniprot_code_handler(key, enr, typ))

        self.db_parser = db_parser
        self.xrefEnum = xrefEnum
        self.ENS_RE = re.compile(r'ENS(?P<species>[A-Z]{0,3})(?P<typ>[GTP])(?P<num>\d{11})')
        self.FB_RE = re.compile(r'FB(?P<typ>[gnptr]{2})(?P<num>\d{7})')
        self.NCBI_RE = re.compile(r'[A-Z]{3}\d{5}\.\d$')
        self.GO_RE = re.compile(r'GO:(?P<termNr>\d{7})')
        self.quote_re = re.compile(r'([[,])([\w_:]+)([,\]])')

    def flush_buffers(self):
        common.package_logger.info('flushing xrefs, go and description buffers')
        if len(self.xrefs) > 0:
            self.xref_tab.append(self.xrefs)
            self.xrefs = []
        if len(self.go) > 0:
            self.go_tab.append(self.go)
            self.go = []

    def _add_to_xrefs(self, eNr, enum_nr, key):
        if not isinstance(eNr, int):
            raise ValueError('eNr is of wrong type:' + str(eNr))
        self.xrefs.append((eNr, enum_nr, key.encode('utf-8'),))
        if len(self.xrefs) > 5e6:
            self.flush_buffers()

    def key_value_handler(self, key, eNr, enum_nr):
        """basic handler that simply adds a key (the xref) under a given enum_nr"""
        self._add_to_xrefs(eNr, enum_nr, key)

    def multi_key_handler(self, multikey, eNr, enum_nr):
        """try to split the myltikey field using '; ' as a delimiter and add each
        part individually under the passed enum_nr id type."""
        for key in multikey.split('; '):
            pos = key.find('.Rep')
            if pos > 0:
                key = key[0:pos]
            self._add_to_xrefs(eNr, enum_nr, key)

    def assign_source_handler(self, multikey, eNr):
        """handler that splits the multikey field at '; ' locations and 
        tries to guess for each part the id_type. If a type could be
        identified, it is added under with this id type, otherwise left out."""
        for key in multikey.split('; '):
            ens_match = self.ENS_RE.match(key)
            if not ens_match is None:
                typ = ens_match.group('typ')
                if typ == 'P':
                    enum_nr = self.xrefEnum['Ensembl Protein']
                elif typ == 'G':
                    enum_nr = self.xrefEnum['Ensembl Gene']
                elif typ == 'T':
                    enum_nr = self.xrefEnum['Ensembl Transcript']
                common.package_logger.debug(
                    'ensembl: ({}, {}, {})'.format(key, typ, enum_nr))
                self._add_to_xrefs(eNr, enum_nr, key)

            for enum, regex in {'FlyBase': self.FB_RE, 'NCBI': self.NCBI_RE}.items():
                match = regex.match(key)
                if not match is None:
                    enum_nr = self.xrefEnum[enum]
                    self._add_to_xrefs(eNr, enum_nr, key)

    def go_handler(self, gos, enr):
        """parse go annotations and add them to the go buffer"""
        for t in gos.split('; '):
            t = t.strip()
            try:
                term, rem = t.split('@')
            except ValueError as e:
                common.package_logger.warning('cannot parse GO annotation: ' + t)
                continue

            term_match = self.GO_RE.match(term)
            termNr = int(term_match.group('termNr'))
            rem = rem.replace('{', '[')
            rem = rem.replace('}', ']')
            rem = self.quote_re.sub('\g<1>"\g<2>"\g<3>', rem)
            for evi, refs in eval(rem):
                for ref in refs:
                    self.go.append((enr, termNr, evi, ref.encode('utf-8')))
            if len(self.go) > 5e6:
                self.flush_buffers()

    def description_handler(self, de, eNr):
        self.desc_manager.add_description(eNr, de)

    def remove_uniprot_code_handler(self, multikey, eNr, enum_nr):
        """remove the species part (sep by '_') of a uniprot long accession to the short acc"""
        common.package_logger.debug(
            'remove_uniprot_code_handler called ({}, {},{})'.format(multikey, eNr, enum_nr))
        for key in multikey.split('; '):
            pos = key.find('_')
            if pos > 0:
                self._add_to_xrefs(eNr, enum_nr, key[0:pos])
            else:
                self._add_to_xrefs(eNr, enum_nr, key)


class IndexedServerParser(object):
    def __init__(self, fh=None):
        """Initializes a Parser for SGML formatted IndexedServer file

        :param fh: file handle of the IndexServer file. If not provided,
                   the file ServerIndexed.db in DARWIN_BROWSERDATA_PATH is
                   used instead.
        """
        self.do_close_at_end = False
        if fh is None:
            fh = open(os.path.join(os.getenv('DARWIN_BROWSERDATA_PATH'), 'ServerIndexed.db'), 'r')
            self.do_close_at_end = True
        self.doc = fh
        self.tag_handlers = collections.defaultdict(list)

    def add_tag_handler(self, tag, handler):
        self.tag_handlers[tag].append(handler)
        common.package_logger.debug('# handlers for {}: {}'.format(tag, len(self.tag_handlers[tag])))

    def parse_entrytags(self):
        """ AC, CHR, DE, E, EMBL, EntrezGene, GI, GO, HGNC_Name, HGNC_Sym, 
        ID, InterPro, LOC, NR , OG, OS, PMP, Refseq_AC, Refseq_ID, SEQ, 
        SwissProt, SwissProt_AC, UniProt/TrEMBL, WikiGene, flybase_transcript_id
        """
        eNr = 0
        for line in self.doc:
            line = line.strip()
            if not line.startswith('<E>'):
                common.package_logger.debug('skipping line:' + line)
                continue

            eNr += 1
            common.package_logger.debug('entry {}: {}'.format(eNr, line))
            entry = lxml.html.fragment_fromstring(line)
            for tag, handlers in self.tag_handlers.items():
                common.package_logger.debug('tag {} ({} handlers)'.format(tag, len(handlers)))
                tag_text = [t.text for t in entry.findall('./' + tag.lower())]
                for value in tag_text:
                    # common.package_logger.debug('value of tag: {}'.format(value.encode('utf-8')))
                    if value is None:
                        continue
                    for handler in handlers:
                        handler(value, eNr)
                        # common.package_logger.debug('called handler {} with ({},{})'.format(
                        #    handler, value.encode('utf-8'), eNr))
        if self.do_close_at_end:
            self.doc.close()


DomainDescription = collections.namedtuple('DomainDescription', tables.dtype_from_descr(tablefmt.DomainDescriptionTable).names)
class CathDomainNameParser(object):
    re_pattern = re.compile(r'(?P<id>[0-9.]*)\s{3,}\w{7}\s{3,}:\s*(?P<desc>.*)')
    source = b'CATH/Gene3D'

    def __init__(self, url):
        self.fname = download_url_if_not_present(url)

    def parse(self):
        open_lib = gzip.open if self.fname.endswith('.gz') else open
        with open_lib(self.fname, 'rt') as fh:
            for line in fh:
                match = self.re_pattern.match(line)
                if match is not None:
                    yield DomainDescription(DomainId=match.group('id').encode('utf-8'),
                                            Source=self.source,
                                            Description=match.group('desc').encode('utf-8'))


class PfamDomainNameParser(CathDomainNameParser):
    re_pattern = re.compile(r'(?P<id>\w*)\t\w*\t\w*\t\w*\t(?P<desc>.*)')
    source = b'Pfam'


def getDebugLogger():
    import logging

    log = logging.getLogger('pyoma')
    log.setLevel(logging.DEBUG)
    logHandler = logging.StreamHandler()
    logHandler.setLevel(logging.DEBUG)
    logHandler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    log.addHandler(logHandler)
    return log


def main(name="OmaServer.h5"):
    log = getDebugLogger()
    x = DarwinExporter(name, logger=log)
    x.add_version()
    x.add_species_data()
    x.add_orthologs()
    x.add_same_species_relations()
    x.add_proteins()
    x.add_hogs()
    x.add_xrefs()
    x.add_domain_info(only_pfam_or_cath_domains(
        iter_domains('ftp://ftp.biochem.ucl.ac.uk/pub/gene3d_data/CURRENT_RELEASE/mdas.csv.gz')))
    x.add_domainname_info(itertools.chain(
        CathDomainNameParser('ftp://ftp.biochem.ucl.ac.uk/pub/cath/latest_release/CathNames').parse(),
        PfamDomainNameParser('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz').parse()))
    x.add_canonical_id()
    x.close()

    x = DarwinExporter(name, logger=log)
    x.create_indexes()
    x.close()

