import random
import unittest
import os
import Bio.Seq
import pyoma.browser.db as pyomadb


class DatabaseChecks(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            path = os.environ['PYOMA_DB2CHECK']
        except KeyError:
            raise unittest.SkipTest("No database specified in PYOMA_DB2CHECK")

        cls.db = pyomadb.Database(path)


    def translated_cdna_match_protein_sequence(self, cdna, prot):
        cdna = cdna.replace('X', 'N')
        for tab in range(1, 22):
            tab_ok = True
            trans = Bio.Seq.translate(cdna, table=tab)
            if not 3 >= len(trans) - len(prot) >= 0:
                return False
            for pos, (trans_aa, prot_aa) in enumerate(zip(trans, prot)):
                if trans_aa == prot_aa or trans_aa == 'X' or prot_aa == 'X':
                    continue
                elif prot_aa == 'M' and pos == 0 and trans_aa != '*':
                    continue
                else:
                    tab_ok = False
                    break
            if tab_ok:
                return True

    def test_cdna_and_protein_sequence_match(self):
        """test translated cdna sequence and protein sequence match.

        This is done for a random sample of 1000 entries"""
        SAMPLES = 1000
        nr_entries = self.db.id_resolver.max_entry_nr
        for entry_nr in random.sample(range(nr_entries+1), SAMPLES):
            with self.subTest(entry_nr=entry_nr):
                cdna = self.db.get_cdna(entry_nr).decode()
                prot = self.db.get_sequence(entry_nr).decode()
                self.assertTrue(self.translated_cdna_match_protein_sequence(cdna, prot))

    def test_increasing_offsets(self):
        entry_tab = self.db.get_hdf5_handle().get_node('/Protein/Entries')
        seq_off = -1
        cds_off = -1
        for row in entry_tab:
            self.assertLess(seq_off, row['SeqBufferOffset'], "SeqBufferOffset decreases in row {}: {} vs {}"
                            .format(row.nrow, seq_off, row['SeqBufferOffset']))
            self.assertLess(cds_off, row['CDNABufferOffset'], "CDNABufferOffset decreases in row {}: {} vs {}"
                            .format(row.nrow, seq_off, row['CDNABufferOffset']))
            seq_off = row['SeqBufferOffset']
            cds_off = row['CDNABufferOffset']
