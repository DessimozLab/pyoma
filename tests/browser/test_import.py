from __future__ import absolute_import
import unittest
import os
import shutil
import tempfile
import json
import numpy
from hashlib import md5

from pyoma.browser.convert import callDarwinExport, DarwinExporter

def store_in_json(data, fn):
    os.mkdir(os.path.dirname(fn))
    with open(fn, 'w') as fd:
        json.dump(data, fd)



class ImportIntegrationBase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.old_env = {(z, os.getenv(z, None)) for z in ('DARWIN_NETWORK_SCRATCH_PATH', 'DARWIN_BROWSERDATA_PATH')}
        os.environ['DARWIN_NETWORK_SCRATCH_PATH'] = cls.tmpdir
        os.environ['DARWIN_BROWSERDATA_PATH'] = '/pub/projects/cbrg-oma-browser/Test.Jul2014/data'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(os.environ['DARWIN_NETWORK_SCRATCH_PATH'])
        for var, val in cls.old_env:
            if val is None:
                os.unsetenv(var)
            else:
                os.environ[var] = val

    def setUp(self):
        fname = tempfile.mktemp(".h5", "OmaServer-", self.tmpdir)
        self.darwin_exporter = DarwinExporter(fname)

    def tearDown(self):
        fn = self.darwin_exporter.h5.filename
        self.darwin_exporter.close()
        os.remove(fn)


class GenomeDirectImport_Test(ImportIntegrationBase):

    def compare_genomes_tab(self, data):
        self.darwin_exporter.add_species_data()
        gstab = self.darwin_exporter.h5.get_node('/Genome')
        self.assertEqual(len(gstab), len(data['GS']), 'unexpected number of genomes')
        for genome in data['GS']:
            gs = gstab.read_where('UniProtSpeciesCode == {}'.format(genome[1].encode('utf-8')))
            for key in ((2, 'TotEntries'), (3, 'TotAA'), (0, 'NCBITaxonId'), (5, 'SciName')):
                expected = genome[key[0]]
                if isinstance(expected, str):
                    expected = expected.encode('utf-8')
                self.assertEqual(gs[key[1]], expected, "data doesn't match for {}: {} vs {}"
                                 .format(key[0], gs[key[1]], expected))
        taxtab = self.darwin_exporter.h5.get_node('/Taxonomy')
        all_taxlevels = taxtab[:]
        self.assertFalse(numpy.where(all_taxlevels['NCBITaxonId'] == all_taxlevels['ParentTaxonId'])[0].any(),
                         'must not have taxlevel pointing to itself')

    def test_load_species_from_json(self):
        data = self.darwin_exporter.call_darwin_export('GetGenomeData();')
        json_fname = os.path.join(self.tmpdir, "pyoma", "gs.json")
        store_in_json(data, json_fname)
        self.compare_genomes_tab(data)
        os.remove(json_fname)

    def test_load_species_from_darwin(self):
        data = self.darwin_exporter.call_darwin_export('GetGenomeData();')
        self.compare_genomes_tab(data)


class ProteinImportViaJson(ImportIntegrationBase):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        data = callDarwinExport('GetGenomeData();')
        json_fname = os.path.join(cls.tmpdir, "pyoma", "gs.json")
        store_in_json(data, json_fname)

    def test_add_proteins(self):
        self.darwin_exporter.add_species_data()
        self.darwin_exporter.add_proteins()

        entry_tab = self.darwin_exporter.h5.get_node('/Protein/Entries')
        sequence_tab = self.darwin_exporter.h5.get_node('/Protein/SequenceBuffer')
        for i, e in enumerate(entry_tab):
            self.assertEqual(i+1, e['EntryNr'], 'entries are not ordered: {} - {}'.format(i, e['EntryNr']))
            seq = sequence_tab[e['SeqBufferOffset']:e['SeqBufferOffset']+e['SeqBufferLength']-1].tostring()
            self.assertEqual(md5(seq).hexdigest(), e['MD5ProteinHash'].decode(),
                             'sequence hashes disagree for {}'.format(e['EntryNr']))

    def test_get_version(self):
        version = self.darwin_exporter.get_version()
        self.assertIn('Test', version)
        self.darwin_exporter.add_version()
        self.assertEqual(version, self.darwin_exporter.h5.get_node_attr('/', 'oma_version'))

    def test_add_orthologs_from_darwin(self):
        pass






