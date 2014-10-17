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
import lxml.html
import collections

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
    EntryNr = tables.UInt32Col(pos=1)
    XRefSource = tables.EnumCol(
        tables.Enum(['UniProtKB/SwissProt', 'UniProtKB/TrEMBL','n/a','EMBL',
            'Ensembl Gene', 'Ensembl Transcript','Ensembl Protein',
            'RefSeq','EntrezGene','GI', 'WikiGene', 'IPI', 'Description',
            'SourceID', 'SourceAC','PMP', 'NCBI', 'FlyBase']),
        'n/a', base='uint8', pos=2)
    XRefId = tables.StringCol(255,pos=3)

class GeneOntologyTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    TermNr = tables.UInt32Col(pos=2)
    Evidence = tables.StringCol(3,pos=3)
    Reference = tables.StringCol(255,pos=4)

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
    fn, off1, off2 = args
    return numpy.genfromtxt(fn, dtype=['i8','i8','S3'], 
            delimiter='\t',
            usecols=(0,1,3),
            converters={0: lambda nr: int(nr)+off1,
                        1: lambda nr: int(nr)+off2})

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
        jobsArgs.append((fn, off1, off2))

    allPairs = []
    pool = mp.Pool()
    allPairs = pool.map(load_tsv_to_numpy, jobsArgs)

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

                vpTab = self.h5.create_table(vpGrp, genome, VPairsTable,
                    expectedrows=len(data))
                enum = vpTab.get_enum('RelType')
                for row in data:
                    row[2] = enum[row[2]]

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
                levels = hog_converter.convert_file(os.path.join(root,fn))
                hogTab.append(levels)
        hog_converter.write_hogs()
        hogTab.flush()
        self.logger.info('createing indexes for HOGs')
        hogTab.cols.Fam.create_index()
        hogTab.cols.ID.create_index()
        entryTab.cols.OmaHOG.create_csindex()

    def add_xrefs(self):
        self.logger.info('start parsing ServerIndexed to extract XRefs and GO annotations')
        db_parser = IndexedServerParser()
        xref_importer = XRefImporter(db_parser)
        db_parser.parse_entrytags()

        self.logger.info('write XRefs / GO to disc')
        xref_tab = self.h5.create_table('/','XRef', XRefTable,
                'Cross-references of proteins to external ids / descriptions',
                expectedrows=len(xref_importer.xrefs))
        self._write_to_table(xref_tab, xref_importer.xrefs)
        go_tab = self.h5.create_table('/','GeneOntology', GeneOntologyTable,
                'Gene Ontology annotations', expectedrows=len(xref_importer.go))
        self._write_to_table(go_tab, xref_importer.go)
        
        self.logger.info('creating index for xrefs (EntryNr and XRefId)')
        xref_tab.cols.EntryNr.create_csindex()
        xref_tab.cols.XRefId.create_csindex()

        self.logger.info('creating index for go (EntryNr and TermNr)')
        go_tab.cols.EntryNr.create_csindex()
        go_tab.cols.TermNr.create_index()

        
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


class XRefImporter(object):
    def __init__(self, db_parser):
        self.xrefs = []
        self.go = []
        xrefEnum = XRefTable.columns.get('XRefSource').enum
        for tag in ['GI','EntrezGene','WikiGene','IPI']:
            db_parser.add_tag_handler(
                    tag, 
                    lambda key, enr, typ=xrefEnum[tag]: self.key_value_handler(key, enr, typ))
        db_parser.add_tag_handler('Refseq_ID', 
            lambda key, enr: self.key_value_handler(key, enr, xrefEnum['RefSeq']))
        db_parser.add_tag_handler('SwissProt', 
            lambda key, enr: self.key_value_handler(key, enr, xrefEnum['UniProtKB/SwissProt']))
        db_parser.add_tag_handler('DE', 
            lambda key, enr: self.key_value_handler(key, enr, xrefEnum['Description']))
        db_parser.add_tag_handler('GO',self.go_handler)
        db_parser.add_tag_handler('ID', self.assign_source_handler)
        db_parser.add_tag_handler('AC', self.assign_source_handler)

        db_parser.add_tag_handler('ID', 
            lambda key, enr: self.multi_key_handler(key, enr, xrefEnum['SourceID']))
        db_parser.add_tag_handler('AC', 
            lambda key, enr: self.multi_key_handler(key, enr, xrefEnum['SourceAC']))
        for tag in ['PMP','EMBL']:
            db_parser.add_tag_handler(
                    tag, 
                    lambda key, enr, typ=xrefEnum[tag]: self.multi_key_handler(key, enr, typ))
        for tag in ['SwissProt_AC','UniProt']: #UniProt/TrEMBL tag is cut to UniProt!
            db_parser.add_tag_handler(tag,
                lambda key, enr, typ=xrefEnum['UniProtKB/TrEMBL']: 
                    self.remove_uniprot_code_handler(key, enr,typ))

        self.db_parser = db_parser
        self.xrefEnum = xrefEnum
        self.ENS_RE=re.compile(r'ENS(?P<species>[A-Z]{0,3})(?P<typ>[GTP])(?P<num>\d{11})')
        self.FB_RE=re.compile(r'FB(?P<typ>[gnptr]{2})(?P<num>\d{7})')
        self.NCBI_RE=re.compile(r'[A-Z]{3}\d{5}\.\d$')
        self.GO_RE=re.compile(r'GO:(?P<termNr>\d{7})')
        self.quote_re=re.compile(r'([[,])([\w_:]+)([,\]])')

    def key_value_handler(self, key, eNr, enum_nr):
        self.xrefs.append((eNr, enum_nr, key,))


    def multi_key_handler(self, multikey, eNr, enum_nr):
        for key in multikey.split('; '):
            pos = key.find('.Rep')
            if pos>0: 
                key=key[0:pos]
            self.xrefs.append((eNr, enum_nr, key,))

    def assign_source_handler(self, multikey, eNr):
        for key in multikey.split('; '):
            ens_match = self.ENS_RE.match(key)
            if not ens_match is None:
                typ = ens_match.group('typ')
                if typ=='P':
                    enum_nr = self.xrefEnum['Ensembl Protein']
                elif typ=='G':
                    enum_nr = self.xrefEnum['Ensembl Gene']
                elif typ=='T':
                    enum_nr = self.xrefEnum['Ensembl Transcript']
                common.package_logger.debug(
                   'ensembl: ({}, {}, {})'.format(key,typ, enum_nr))
                self.xrefs.append((eNr,enum_nr,key))
            
            for enum, regex in {'FlyBase':self.FB_RE, 'NCBI':self.NCBI_RE}.items():
                match = regex.match(key)
                if not match is None:
                    enum_nr = self.xrefEnum[enum]
                    self.xrefs.append((eNr, enum_nr, key))


    def go_handler(self, gos, enr):
        for t in gos.split('; '):
            term, rem = t.split('@')
            term_match = self.GO_RE.match(term)
            termNr = int(term_match.group('termNr'))
            rem = rem.replace('{','[')
            rem = rem.replace('}',']')
            rem = self.quote_re.sub('\g<1>"\g<2>"\g<3>',rem)
            for evi, refs in eval(rem):
                for ref in refs:
                    self.go.append((enr, termNr, evi, ref))


    def remove_uniprot_code_handler(self, multikey, eNr, enum_nr):
        common.package_logger.debug(
                'remove_uniprot_code_handler called ({}, {},{})'.format(multikey,eNr, enum_nr))
        for key in multikey.split('; '):
            pos = key.find('_')
            if pos>0:
                self.xrefs.append((eNr, enum_nr, key[0:pos]))
            else:
                self.xrefs.append((eNr, enum_nr, key))


class IndexedServerParser(object):
    def __init__(self, fn=None):
        if fn is None:
            fn = os.path.join(os.getenv('DARWIN_BROWSERDATA_PATH'), 'ServerIndexed.db')
        self.elems = lxml.html.parse(fn).getroot()
        self.tag_handlers = collections.defaultdict(list)

    def add_tag_handler(self, tag, handler):
        self.tag_handlers[tag].append(handler)
        common.package_logger.debug('# handlers for {}: {}'.format(tag, len(self.tag_handlers[tag])))
        

    def parse_entrytags(self):
        """ AC, CHR, DE, E, EMBL, EntrezGene, GI, GO, HGNC_Name, HGNC_Sym, 
        ID, InterPro, LOC, NR , OG, OS, PMP, Refseq_AC, Refseq_ID, SEQ, 
        SwissProt, SwissProt_AC, UniProt/TrEMBL, WikiGene, flybase_transcript_id
        """
        for eNr_minus_1, entry in enumerate(self.elems.findall('.//e')):
            eNr = eNr_minus_1+1
            common.package_logger.debug('handling entry {} looking {}'.format(eNr, entry.text))
            for tag, handlers in self.tag_handlers.items():
                common.package_logger.debug('tag {} ({} handlers)'.format(tag, len(handlers)))
                tag_text = [t.text for t in entry.findall('./'+tag.lower())]
                for value in tag_text:
                    common.package_logger.debug('value of tag: {}'.format(value))
                    for handler in handlers:
                       handler(value, eNr)
                       common.package_logger.debug('called handler {} with ({},{})'.format(
                           str(handler), str(value), eNr))




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

