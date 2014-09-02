from __future__ import division
from future.builtins import str
from future.builtins import chr
from future.builtins import range
from future.builtins import object
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
from . import common
from . import locus_parser


# definitions of the pytables formats

class HOGsTable(tables.IsDescription):
    Fam = tables.Int32Col(pos=1)
    ID = tables.StringCol(255,pos=2)
    Level = tables.StringCol(255,pos=3)

class ProteinTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    SeqBufferOffset = tables.UInt32Col(pos=2)
    SeqBufferLength = tables.UInt32Col(pos=3)
    OmaGroup = tables.UInt32Col(pos=4,dflt=0)
    OmaHOG = tables.StringCol(255,pos=5,dflt='')
    Chromosome = tables.StringCol(255,pos=6)
    LocusStart = tables.UInt32Col(pos=7)
    LocusEnd = tables.UInt32Col(pos=8)
    LocusStrand = tables.Int8Col(pos=9,dflt=1)
    AltSpliceVariant = tables.Int32Col(pos=10,dflt=0)

class LocusTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    Start = tables.UInt32Col(pos=2)
    End = tables.UInt32Col(pos=3)
    Strand = tables.Int8Col(pos=4)

class VPairsTable(tables.IsDescription):
    EntryNr1 = tables.UInt32Col(pos=0)
    EntryNr2 = tables.UInt32Col(pos=1)
    RelType = tables.EnumCol(
        tables.Enum(['1:1','1:n','m:1','m:n','n/a']),
        'n/a', base='uint8', pos=2)

class XRefTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=0)
    XRefSource = tables.EnumCol(
        tables.Enum(['UniProtKB/SwissProt', 'UniProtKB/TrEMBL','n/a','EMBL',
            'Ensembl Gene', 'Ensembl Transcript','Ensembl Protein',
            'RefSeq','EntrezGene','GI', 'WikiGene', 'IPI', 'Description',
            'SourceID', 'SourceAC']),
        'n/a', base='uint8', pos=1)
    XRefId = tables.StringCol(30,pos=2)

class GenomeTable(tables.IsDescription):
    NCBITaxonId = tables.UInt32Col(pos=0)
    UniProtSpeciesCode = tables.StringCol(5,pos=1)
    TotEntries = tables.UInt32Col(pos=2)
    TotAA = tables.UInt32Col(pos=3)
    EntryOff = tables.UInt32Col(pos=4)
    SciName = tables.StringCol(255,pos=5)
    Release = tables.StringCol(255,pos=6)

class TaxonomyTable(tables.IsDescription):
    NCBITaxonId = tables.UInt32Col(pos=0)
    ParentTaxonId = tables.UInt32Col(pos=1)
    Name = tables.StringCol(255,pos=2)


class DarwinException(Exception):
    pass


def callDarwinExport(func):
    """Function starts a darwin session, loads convert.drw file
    and calls the darwin function passed as argument. The output
    is expected to be written by darwin in json format into the
    file specified by 'outfn'.
    This function returns the parsed json datastructure"""

    tmpfile = "/tmp/darwinExporter_%d.dat"%(os.getpid())
    drwCodeFn = os.path.abspath(
        os.path.splitext(__file__)[0]+'.drw')
    try:
        with open(os.devnull,'w') as DEVNULL:
            p = subprocess.Popen(['darwin','-q','-E','-B'], stdin=subprocess.PIPE,
                stderr=subprocess.PIPE, stdout=DEVNULL)
            drwCmd = "outfn := '{}': ReadProgram('{}'): {}; done;".format(
                tmpfile, drwCodeFn, func).encode('utf-8')
            p.communicate(input=drwCmd)
            if p.returncode>0:
                raise DarwinException( p.stderr.read() )

        trans_tab = "".join(str(chr(x)) for x in range(128))+" "*128
        with open(tmpfile,'r') as jsonData:
            rawdata = jsonData.read()
            data = json.loads(rawdata.translate(trans_tab));
    finally:
        silentremove(tmpfile)
    return data


def silentremove(filename):
    """Function to remove a given file. No exception is raised if the
    file does not exist. Other errors are passed to the user."""
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured

def load_tsv_to_numpy(args):
    fn, off1, off2, swap = args
    relEnum = VPairsTable.columns['RelType'].enum._names
    relEnum['n:1'] = relEnum['m:1']
    relEnum['1:m'] = relEnum['1:n']
    relEnum['n:m'] = relEnum['m:n']
    read_dir = -1 if swap else 1
    dtype = [('EntryNr1','i4'),('EntryNr2','i4'),('Score','i4'),('RelType','i1')]
    for curNr, curFn in enumerate([fn, fn.replace('.ext.','.')]):
        try:
            augData = numpy.genfromtxt(curFn, dtype=dtype,
                names=[_[0] for _ in dtype],
                delimiter='\t',
                usecols=(0,1,2,3),
                converters={'EntryNr1': lambda nr: int(nr)+off1,
                        'EntryNr2': lambda nr: int(nr)+off2,
                        'RelType': lambda rel: relEnum[rel[::read_dir].decode()]})
            break;
        except OSError as e:
            if curNr<1:
                pass
            else:
                raise e

    ret_cols = ['EntryNr1','EntryNr2','RelType']
    if swap:
        reversed_cols = tuple(augData.dtype.names[z] for z in (1,0,2,3))
        augData.dtype.names = reversed_cols
    return augData[ret_cols]

def read_vps_from_tsv(gs, ref_genome):
    ref_genome_idx = gs.get_where_list('(UniProtSpeciesCode=={})'.
            format(ref_genome))[0]
    jobsArgs = [] 
    for g in range(len(gs)):
        if g==ref_genome_idx:
            continue
        g1, g2 = sorted((g, ref_genome_idx,))
        off1, off2 = gs.read_coordinates(numpy.array((g1,g2)), 'EntryOff')
        fn = os.path.join(os.environ['DARWIN_OMADATA_PATH'],'Phase4',
                gs.cols.UniProtSpeciesCode[g1].decode(), 
                gs.cols.UniProtSpeciesCode[g2].decode()+".ext.orth.txt.gz")
        jobsArgs.append((fn, off1, off2, g1!=ref_genome_idx))

    pool = mp.Pool()
    allPairs = pool.map(load_tsv_to_numpy, jobsArgs)
    pool.close()

    return numpy.lib.recfunctions.stack_arrays(allPairs, usemask=False) 
    


class DataImportError(Exception):
    pass


class DarwinExporter(object):
    def __init__(self, path, logger=None):
        self.logger = logger if logger is not None else common.package_logger
        fn = os.path.normpath(os.path.join(
            os.environ['DARWIN_BROWSERDATA_PATH'],
            path))
        mode = 'append' if os.path.exists(fn) else 'write'
        self._compr = tables.Filters(complevel=6,complib='zlib', fletcher32=True)
        self.h5 = tables.open_file(fn, mode=mode[0], filters=self._compr)
        self.logger.info("opened %s in %s mode, options %s"%(
            fn, mode, str(self._compr)))
        self.h5.root._f_setattr('convertion_start', time.strftime("%c"))

    def add_version(self):
        # FIXME: proper implementation of version
        version = 'Sep 2014'
        self.h5.root._f_setattr('oma_version', version)

    def add_species_data(self):
        cache_file = os.path.join(
            os.environ['DARWIN_NETWORK_SCRATCH_PATH'],
            'pyoma','gs.json')
        if os.path.exists(cache_file):
            with open(cache_file, 'r') as fd:
                data = json.load(fd)
        else:
            data = callDarwinExport('GetGenomeData();')
        gstab = self.h5.create_table('/', 'Genome', GenomeTable, 
            expectedrows=len(data['GS']))
        self._write_to_table(gstab, data['GS'])
        gstab.cols.NCBITaxonId.create_csindex(filters=self._compr)
        gstab.cols.UniProtSpeciesCode.create_csindex(filters=self._compr)
        gstab.cols.EntryOff.create_csindex(filters=self._compr)

        taxtab = self.h5.create_table('/', 'Taxonomy', TaxonomyTable, 
            expectedrows=len(data['Tax']))
        self._write_to_table(taxtab, data['Tax'])
        taxtab.cols.NCBITaxonId.create_csindex(filters=self._compr)

    def add_orthologs(self):
        try:
            vpGrp = self.h5.get_node('/VPairs')
        except tables.NoSuchNodeError:
            vpGrp = self.h5.create_group('/','VPairs',title='Pairwise Orthologs')

        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            if not '/'+genome in self.h5:
                cache_file = os.path.join(
                    os.environ['DARWIN_NETWORK_SCRATCH_PATH'],
                    'pyoma','vps','{}.json'.format(genome))
                if os.path.exists(cache_file):
                    with open(cache_file, 'r') as fd:
                        data = json.load(fd)
                elif not os.getenv('DARWIN_OMADATA_PATH') is None:
                    # try to read from Phase4 in parallel.
                    data = read_vps_from_tsv(self.h5.root.Genome, 
                                             genome.encode('utf-8'))
                else:
                    # fallback to read from VPsDB
                    data = callDarwinExport('GetVPsForGenome({})'.
                            format(genome))

                self.logger.debug(data[1])
                vpTab = self.h5.create_table(vpGrp, genome, VPairsTable,
                    expectedrows=len(data))
                if isinstance(data, list):
                    enum= vpTab.get_enum('RelType')
                    for row in data:
                        row[2] = enum[row[2].decode()]

                self._write_to_table(vpTab, data)


    def add_proteins(self):
        gsNode = self.h5.get_node('/Genome')
        nrProt = sum(gsNode.cols.TotEntries)
        nrAA = sum(gsNode.cols.TotAA)
        try:
            protGrp = self.h5.get_node('/Protein')
        except tables.NoSuchNodeError:
            protGrp = self.h5.create_group('/','Protein')
        protTab = self.h5.create_table(protGrp, 'Entries', ProteinTable,
            expectedrows=nrProt)
        seqArr = self.h5.create_earray(protGrp,'SequenceBuffer', 
            tables.StringAtom(1),(0,),'concatenated protein sequences',
            expectedrows=nrAA+nrProt)

        seqOff = 0
        for gs in gsNode.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            cache_file = os.path.join(
                os.environ['DARWIN_NETWORK_SCRATCH_PATH'],
                'pyoma','prots','{}.json'.format(genome))
            if os.path.exists(cache_file):
                with open(cache_file,'r') as fd:
                    data = json.load(fd)
            else:
                data = callDarwinExport('GetProteinsForGenome(%s)'%(genome))

            if len(data['seqs']) != gs['TotEntries']:
                raise DataImportError('number of entries (%d) does '
                    'not match number of seqs (%d) for %s'%
                    (len(data['seqs']), gs['TotEntries'], genome))

            locTab = self.h5.create_table('/Protein/Locus', 
                genome, LocusTable, createparents=True, 
                expectedrows=gs['TotEntries']*4)

            for nr in range(gs['TotEntries']):
                eNr = data['off']+nr+1
                protTab.row['EntryNr'] = eNr 
                protTab.row['OmaGroup'] = data['ogs'][nr]
                # add ' ' after each sequence (Ascii is smaller than 
                # any AA, allows to build PAT array with split between 
                # sequences.
                seqLen = len(data['seqs'][nr])+1
                protTab.row['SeqBufferOffset'] = seqOff
                protTab.row['SeqBufferLength'] = seqLen
                seqNumpyObj = numpy.ndarray((seqLen,), 
                    buffer=(data['seqs'][nr]+" ").encode('utf-8'),
                    dtype=tables.StringAtom(1))
                seqArr.append(seqNumpyObj)
                seqOff += seqLen
                protTab.row['Chromosome'] = data['chrs'][nr]
                protTab.row['AltSpliceVariant'] = data['alts'][nr]

                locus_str = data['locs'][nr]
                try:
                    locus_tab = locus_parser.parse(eNr, locus_str)
                    locTab.append(locus_tab)

                    protTab.row['LocusStart'] = locus_tab.Start.min()
                    protTab.row['LocusEnd'] = locus_tab.End.max()
                    protTab.row['LocusStrand'] = locus_tab[0].Strand
                except ValueError as e:
                    self.logger.warning(e)
                protTab.row.append()
            protTab.flush()
            seqArr.flush()
            for n in (protTab, seqArr, locTab):
                self.logger.info('worte %s: compression ratio %3f%%'%
                    (n._v_pathname, 100*n.size_on_disk/n.size_in_memory))
        protTab.cols.EntryNr.create_csindex(filters=self._compr)


    def _write_to_table(self, tab, data):
        tab.append(data)
        self.logger.info('wrote %s : compression ratio %.3f%%'%
                (tab._v_pathname, 100*tab.size_on_disk/tab.size_in_memory))

    def add_hogs(self):
        hog_path = os.path.normpath(os.path.join(
            os.environ['DARWIN_BROWSERDATA_PATH'],
            '..','downloads','HOGs'))
        entryTab = self.h5.get_node('/Protein/Entries')
        tree_filename = os.path.join(
                os.environ['DARWIN_BROWSERDATA_PATH'],
                'speciestree.nwk')
        hog_converter = HogConverter(entryTab)
        hog_converter.attach_newick_taxonomy(tree_filename)
        hogTab = self.h5.create_table('/', 'HogLevel', HOGsTable, 
            'nesting structure for each HOG', expectedrows=1e8)
        for root, dirs, filenames in os.walk(hog_path):
            for fn in filenames:
                try:
                    levels = hog_converter.convert_file(os.path.join(root,fn))
                    hogTab.append(levels)
                except Exception as e:
                    self.logger.error('an error occured while processing '+fn+':')
                    self.logger.exception(e)
                    
        hog_converter.write_hogs()
        hogTab.flush()
        self.logger.info('createing indexes for HOGs')
        hogTab.cols.Fam.create_index()
        hogTab.cols.ID.create_index()
        entryTab.cols.OmaHOG.create_csindex()
        
    def close(self):
        self.h5.root._f_setattr('conversion_end', time.strftime("c"))
        self.h5.close()
        self.logger.info('closed %s'%(self.h5.filename))


class GroupAnnotatorInclGeneRefs(familyanalyzer.GroupAnnotator):
    def _annotateGroupR(self, node, og, idx=0):
        if familyanalyzer.OrthoXMLQuery.is_geneRef_node(node):
            node.set('og', og)
        else:
            super()._annotateGroupR(node, og, idx)


class HogConverter(object):
    def __init__(self, entry_tab):
        self.fam_re = re.compile(r'HOG:(?P<fam_nr>\d+)')
        self.hogs = numpy.chararray(shape=(len(entry_tab)+1,), itemsize=255)
        self.hogs[:] = ''
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
            levs.extend([(fam_nr, n.getparent().get('og'),n.get('value'),) 
                for n in p._findSubNodes('property')
                    if n.get('name')=="TaxRange"])

        geneNodes = p.root.findall('.//{{{ns0}}}geneRef'.
              format(**familyanalyzer.OrthoXMLParser.ns))
        for x in geneNodes:
            self.hogs[int(x.get('id'))] = x.get('og')

        return levs

    def write_hogs(self):
        self.entry_tab.modify_column(0,len(self.entry_tab),1, self.hogs[1:], 'OmaHOG') 
        self.entry_tab.flush()



def getDebugLogger():
    import logging
    log = logging.getLogger('pyoma')
    log.setLevel(logging.DEBUG)
    logHandler=logging.StreamHandler()
    logHandler.setLevel(logging.DEBUG)
    logHandler.setFormatter(logging.Formatter(
       '%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    log.addHandler(logHandler)
    return log


def main(name="OmaServer.h5"):
    log = getDebugLogger()
    x=DarwinExporter(name, logger=log)
    x.add_version()
    x.add_species_data()
    x.add_orthologs()
    x.add_proteins()
    x.add_hogs()
    x.close()

