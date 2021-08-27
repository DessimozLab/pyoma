import random
import unittest
import os
import Bio.Seq
import Bio.Data.CodonTable
import pyoma.browser.db as pyomadb
import tables
import numpy


class DatabaseChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            path = os.environ["PYOMA_DB2CHECK"]
        except KeyError:
            raise unittest.SkipTest("No database specified in PYOMA_DB2CHECK")

        cls.db = pyomadb.Database(path)

    def translated_cdna_match_protein_sequence(self, cdna, prot):
        cdna = cdna.replace("X", "N")
        for tab in Bio.Data.CodonTable.generic_by_id.keys():
            tab_ok = True
            trans = Bio.Seq.translate(cdna, table=tab)
            if not 3 >= len(trans) - len(prot) >= 0:
                return False
            for pos, (trans_aa, prot_aa) in enumerate(zip(trans, prot)):
                if trans_aa == prot_aa or trans_aa == "X" or prot_aa == "X":
                    continue
                elif prot_aa == "M" and pos == 0 and trans_aa != "*":
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
        for entry_nr in random.sample(range(nr_entries + 1), SAMPLES):
            with self.subTest(entry_nr=entry_nr):
                cdna = self.db.get_cdna(entry_nr).decode()
                prot = self.db.get_sequence(entry_nr).decode()
                self.assertTrue(
                    abs(3 * len(prot) - len(cdna)) < 9,
                    "CDS and protein sequence are too different in len: 3*{}={} vs {}".format(
                        len(prot), 3 * len(prot), len(cdna)
                    ),
                )
                self.assertTrue(self.translated_cdna_match_protein_sequence(cdna, prot))

    def test_increasing_offsets(self):
        entry_tab = self.db.get_hdf5_handle().get_node("/Protein/Entries")
        seq_off = -1
        cds_off = -1
        for row in entry_tab:
            self.assertLess(
                seq_off,
                row["SeqBufferOffset"],
                "SeqBufferOffset decreases in row {}: {} vs {}".format(
                    row.nrow, seq_off, row["SeqBufferOffset"]
                ),
            )
            self.assertLess(
                cds_off,
                row["CDNABufferOffset"],
                "CDNABufferOffset decreases in row {}: {} vs {}".format(
                    row.nrow, seq_off, row["CDNABufferOffset"]
                ),
            )
            seq_off = row["SeqBufferOffset"]
            cds_off = row["CDNABufferOffset"]

    def test_homeology_flag(self):
        genome_tab = self.db.get_hdf5_handle().get_node("/Genome")
        for g in (b"WHEAT", b"GOSHI", b"BRANA"):
            for row in genome_tab.read_where("UniProtSpeciesCode == g"):
                self.assertTrue(
                    row["IsPolyploid"],
                    "{} is not recorded as polyploid genome".format(g),
                )
        for g in (b"YEAST", b"HUMAN", b"PLAF7", b"ARATH", b"MOUSE"):
            for row in genome_tab.read_where("UniProtSpeciesCode == g"):
                self.assertFalse(
                    row["IsPolyploid"],
                    "{} is recorded to be a ployploid genome".format(g),
                )

    def test_synteny_scores_exist(self):
        for g in ("WHEAT", "BRANA", "GOSHI"):
            try:
                t = self.db.get_hdf5_handle().get_node(
                    "/PairwiseRelation/{}/within".format(g)
                )
            except tables.NoSuchNodeError:
                # if species does not exist, we skip - not all datasets will have these genomes
                continue
            syn_col = t.col("SyntenyConservationLocal")
            computed_pairs = numpy.where(syn_col >= 0)
            self.assertLess(
                0, len(computed_pairs[0]), "No synteny values computed for {}".format(g)
            )

    def test_alignment_of_pairwise_orthologs(self):
        SAMPLES = 500
        nr_entries = self.db.id_resolver.max_entry_nr
        seq_search = self.db.seq_search
        for entry_nr1 in random.sample(range(1, nr_entries + 1), 10 * SAMPLES):
            vps = self.db.get_vpairs(entry_nr1)
            if len(vps) == 0:
                continue
            pair = random.choice(vps)
            with self.subTest(pair=pair):
                s1 = self.db.get_sequence(pair["EntryNr1"])
                matches = [(pair["EntryNr2"],)]
                alignment = list(
                    seq_search._align_entries(s1, matches, compute_distance=True)
                )[0]
                score = alignment[1][0]
                dist = alignment[1][2]
                if pair["Score"] < 1:
                    self.assertGreaterEqual(score, 120)
                else:
                    self.assertAlmostEqual(pair["Score"], score, places=1)
                    self.assertAlmostEqual(pair["Distance"], dist, places=1)
            SAMPLES -= 1
            if SAMPLES <= 0:
                break

    def test_cached_orthologs_count(self):
        SAMPLES = 80
        nr_entries = self.db.id_resolver.max_entry_nr
        for entry_nr in random.sample(range(nr_entries + 1), SAMPLES):
            with self.subTest(entry_nr=entry_nr):
                cached_ortholog_cnts = self.db.nr_ortholog_relations(entry_nr)
                vps = self.db.get_vpairs(entry_nr)
                self.assertEqual(len(vps), cached_ortholog_cnts["NrPairwiseOrthologs"])

                entry = self.db.entry_by_entry_nr(entry_nr)
                if entry["OmaGroup"] > 0:
                    oma_grp_memb = self.db.oma_group_members(entry["OmaGroup"])
                    self.assertEqual(
                        len(oma_grp_memb) - 1,
                        cached_ortholog_cnts["NrOMAGroupOrthologs"],
                    )

                if entry["OmaHOG"] != b"":
                    hog_ind_orth = self.db.get_hog_induced_pairwise_orthologs(entry)
                    self.assertEqual(
                        len(hog_ind_orth),
                        cached_ortholog_cnts["NrHogInducedPWOrthologs"],
                    )
