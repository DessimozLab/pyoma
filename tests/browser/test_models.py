from pyoma.browser import models, db
from future.utils import with_metaclass
import sys
import os
import unittest


class ProteinEntryTests(unittest.TestCase):
    db = None

    @classmethod
    def setUpClass(cls):
        fn = os.path.join(os.path.dirname(__file__), 'TestDb.h5')
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
        self.assertLessEqual(protein_entry.ec_content, 1)
        self.assertGreaterEqual(protein_entry.ec_content, 0)


class SingletonTests(unittest.TestCase):

    def test_singleton(self):
        class Foo(with_metaclass(models.Singleton, object)):
            def __init__(self):
                pass

        a = Foo()
        b = Foo()
        self.assertEqual(a, b)