from __future__ import unicode_literals, division, absolute_import
from builtins import bytes, range, str

import numpy
from pyoma.browser import models, db, tablefmt
from future.utils import with_metaclass
import sys
import os
import unittest
from .test_db import find_path_to_test_db


class ProteinEntryTests(unittest.TestCase):
    db = None

    @classmethod
    def setUpClass(cls):
        fn = find_path_to_test_db('TestDb.h5')
        cls.db = db.Database(fn)

    @classmethod
    def tearDownClass(cls):
        cls.db.get_hdf5_handle().close()

    def test_init_from_enr(self):
        protein_entry = models.ProteinEntry.from_entry_nr(self.db, 12)
        self.assertEqual(protein_entry.entry_nr, 12)
        try:
            protein_entry.canonicalid
        except AttributeError:
            self.assertTrue(False, 'entry not properly loaded')

    def test_lazy_eval_of_properties(self):
        protein_entry = models.ProteinEntry.from_entry_nr(self.db, 12)
        self.assertNotIn('sequence', protein_entry.__dict__)
        seq = protein_entry.sequence
        self.assertIn('sequence', protein_entry.__dict__)

    def test_repr_of_ProteinEntry(self):
        protein_entry = models.ProteinEntry.from_entry_nr(self.db, 12)
        if sys.version_info[0] >= 3 and sys.version_info[1] > 2:
            self.assertRegex("{!r}".format(protein_entry), r"<ProteinEntry\(12,.*\)")
        else:
            self.assertRegexpMatches("{!r}".format(protein_entry), r"<ProteinEntry\(12,.*\)")

    def test_len_of_entry(self):
        protein_entry = models.ProteinEntry.from_entry_nr(self.db, 12)
        self.assertEqual(len(protein_entry), len(protein_entry.sequence))

    def test_ec_content(self):
        protein_entry = models.ProteinEntry.from_entry_nr(self.db, 12)
        self.assertLessEqual(protein_entry.gc_content, 1)
        self.assertGreaterEqual(protein_entry.gc_content, 0)

    def test_example_ec_content(self):
        protein_entry = models.ProteinEntry.from_entry_nr(self.db, 12)
        protein_entry.cdna = "GCGAATAT"
        self.assertAlmostEqual(protein_entry.gc_content, 3.0 / 8.0)


class ExonStructureTest(unittest.TestCase):
    def get_exons(self):
        loc_dtype = tablefmt.tables.dtype_from_descr(tablefmt.LocusTable)
        return numpy.array([(1, 500, 510, 1), (1, 600, 610, 1)],
                           dtype=loc_dtype)

    def test_exon_struct_len(self):
        ex = models.ExonStructure(None, self.get_exons())
        self.assertEqual(2, len(ex))

    def test_str_repr(self):
        ex = models.ExonStructure(None, self.get_exons())
        self.assertEqual("join(500..510, 600..610)", str(ex))

    def test_str_repr_if_reverse_complement(self):
        ex_dat = self.get_exons()
        ex_dat['Strand'] = -1
        ex = models.ExonStructure(None, ex_dat)
        self.assertEqual("join(complement(600..610), complement(500..510))", str(ex))


class GeneOntologyAnnotationTests(unittest.TestCase):
    db = None

    @classmethod
    def setUpClass(cls):
        fn = find_path_to_test_db('TestDb.h5')
        cls.db = db.Database(fn)

    @classmethod
    def tearDownClass(cls):
        cls.db.get_hdf5_handle().close()

    def test_from_annotation_table(self):
        annos = self.db.get_gene_ontology_annotations(12)
        goa = models.GeneOntologyAnnotation(self.db, annos[0])
        self.assertEqual(12, goa.entry_nr)
        self.assertIn(goa.aspect, ['molecular_function', 'biological_process','cellular_component'])


class SingletonTests(unittest.TestCase):

    def test_singleton(self):
        class Foo(with_metaclass(models.Singleton, object)):
            def __init__(self):
                pass

        a = Foo()
        b = Foo()
        self.assertEqual(a, b)