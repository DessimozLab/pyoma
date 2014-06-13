import tables
import numpy
import os
import subprocess
import errno
import json
import time
from . import common


# definitions of the pytables formats

class HOGsTable(tables.IsDescription):
    ID = tables.StringCol(255)
    Level = tables.StringCol(255)

class ProteinTable(tables.IsDescription):
    NCBITaxonId = tables.UInt32Col(pos=0)
    EntryNr = tables.UInt32Col(pos=1)
    SeqBufferOffset = tables.UInt32Col(pos=2)
    SeqBufferLength = tables.UInt32Col(pos=3)
    OmaGroup = tables.UInt32Col(pos=4,dflt=0)
    OmaHOG = tables.StringCol(255,pos=5)
    Chromosome = tables.StringCol(255,pos=6)
    Locus = tables.StringCol(255,pos=7)

class VPairsTable(tables.IsDescription):
    EntryNr1 = tables.UInt32Col(pos=0)
    EntryNr2 = tables.UInt32Col(pos=1)
    RelType = tables.EnumCol(
        tables.Enum(['1:1','1:m','n:1','n:m']),
        '1:1', base='uint8', pos=2)

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
    SciName = tables.StringCol(255,pos=3)
    Release = tables.StringCol(255,pos=4)

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
            p.communicate(input="outfn := '%s': ReadProgram('%s'): %s; done;"
                %(tmpfile,drwCodeFn,func));
            if p.returncode>0:
                raise DarwinException( p.stderr.read() )

        trans_tab = "".join(str(unichr(x)) for x in xrange(128))+" "*128
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


class DarwinExporter(object):
    def __init__(self, path, logger=None):
        self.logger = logger if logger is not None else common.package_logger
        fn = os.path.normpath(os.path.join(
            os.environ['DARWIN_BROWSERDATA_PATH'],
            path))
        self._compr = tables.Filters(complevel=6,complib='zlib', fletcher32=True)
        self.h5 = tables.open_file(fn, mode='w', filters=self._compr)
        self.logger.info("start writing to %s, options %s"%(fn, str(self._compr)))
        self.h5.root._f_setattr('convertion_start', time.strftime("%c"))

    def add_version(self):
        version = 'Mar 2014'
        self.h5.root._f_setattr('oma_version', version)

    def add_species_data(self):
        data = callDarwinExport('GetGenomeData')
        gstab = self.h5.create_table('/Genome', GenomeTable, expectedrows=len(data['gs']))
        self._write_to_table(gstab, data['gs'])
        gstab.cols.NCBITaxonId.create_csindex(filters=self._compr)
        gstab.cols.UniProtSpeciesCode.create_csindex(filters=self._compr)

        taxtab = self.h5.create_table('/Taxonomy', TaxonomyTable, expectedrows=len(data['lin']))
        self._write_to_table(taxtab, data['lin'])
        taxtab.cols.NCBITaxonId.create_csindex(filters=self._compr)

    def _write_to_table(self, tab, data):
        tab.append(data)
        self.logger.info('wrote %s : compression ratio %.3f%%'%
                (tab._v_pathname, 100*tab.size_on_disk/tab.size_in_memory))

    def close(self):
        self.h5.root._f_setattr('conversion_end', time.strftime("c"))
        self.h5.close()
        self.logger.info('closed %s'%(self.h5.filename))



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
    x.close()

