from __future__ import absolute_import, unicode_literals, print_function, division

import gzip
import unittest
import os
import shutil
import tempfile
import json
import numpy
import csv
from hashlib import md5
import pyoma.browser.tablefmt as tablefmt
from pyoma.browser.convert import callDarwinExport, DarwinExporter, compute_ortholog_types,\
    load_tsv_to_numpy


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
        test_data_available = False
        for folder in ('/pub/projects/cbrg-oma-browser', '/cs/research/biosciences/oma/oma-server',
                       '/scratch/ul/projects/cdessimo/oma-browser'):
            if os.path.isdir(folder):
                test_data_available = True
                break
        if not test_data_available:
            raise unittest.SkipTest('data not available')
        os.environ['DARWIN_BROWSERDATA_PATH'] = os.path.join(folder, 'Test.Jul2014', 'data')

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


class GenomeDirectImportTest(ImportIntegrationBase):
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
        super(ProteinImportViaJson, cls).setUpClass()
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
        self.assertEqual(version, self.darwin_exporter.h5.get_node_attr('/','oma_version'))

    def test_add_orthologs_from_darwin(self):
        pass


class OrthologyTypeTester(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.mkdtemp()
        self.exp = DarwinExporter(os.path.join(self._tmpdir, 'test.h5'))
        self.pw = self.exp.h5.create_table('/', 'VPairs',
                                           tablefmt.PairwiseRelationTable)

    def test_convert_to_numpytable(self):
        rels = [[1, 5], (2, 6)]
        tab = self.exp._convert_to_numpyarray(rels, self.pw)
        self.assertEqual(len(tab), 2)
        # test that default value are used
        self.assertTrue(numpy.all(tab['Distance'] == -1))

    def test_orthology_type_inference(self):
        genome_offs = numpy.array([0, 4, 9, 19])
        rels = [[1, 5], [1, 6], [1, 10], [1, 20],
                [2, 5], [2, 6], [2, 10], [2, 11],
                [3, 12], [3, 21], [3, 22]]
        expected_rel_type = numpy.array([3, 3, 2, 0, 3, 3, 3, 1, 0, 1, 1])

        reltab = self.exp._convert_to_numpyarray(rels, self.pw)
        compute_ortholog_types(reltab, genome_offs)
        self.assertTrue(numpy.array_equal(reltab['RelType'], expected_rel_type))

    def tearDown(self):
        self.exp.close()
        shutil.rmtree(self._tmpdir)


class TSVOrthologyFileTester(unittest.TestCase):
    def setUp(self):
        with tempfile.NamedTemporaryFile(suffix=".gz", delete=False) as fh:
            self.datafile = fh.name

    def tearDown(self):
        try:
            os.remove(self.datafile)
        except OSError:
            pass

    def put_test_data(self, data):
        with gzip.open(self.datafile, mode='wt', compresslevel=6) as cmp:
            writer = csv.writer(cmp, delimiter='\t')
            for row in data:
                writer.writerow(row)

    def test_regular_example(self):
        data = [(1, 1, 12421, '1:1', .99, 1.523),
                (2, 5, 21214, 'n:1', .94, 12.14),
                (3, 5, 23552, 'n:1', .91, 21.2)]
        self.put_test_data(data)
        res = load_tsv_to_numpy((self.datafile, 0, 100, False))
        self.assertEqual(3, res.size)
        self.assertTrue((numpy.arange(1, 4) == res['EntryNr1']).all())
        self.assertTrue((numpy.array((101, 105, 105)) == res['EntryNr2']).all())
        self.assertTrue((numpy.array((0, 2, 2)) == res['RelType']).all())
        self.assertTrue(((res['Distance'] > 0) & (res['Distance'] < 22)).all())
        self.assertTrue(((res['AlignmentOverlap'] > 0.9) & (res['AlignmentOverlap'] <= 1)).all())

    def test_empty_file(self):
        self.put_test_data([])
        res = load_tsv_to_numpy((self.datafile, 0, 0, False))
        self.assertEqual(0, res.size)

