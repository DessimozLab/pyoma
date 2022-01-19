import unittest
import logging
import time
from unittest.mock import patch

from .test_db import TestWithDbInstance
from pyoma.browser.models import ProteinEntry, HOG
from pyoma.browser.search import (
    OmaGroupSearch,
    HogIDSearch,
    GOSearch,
    TaxSearch,
    SequenceSearch,
    ECSearch,
    XRefSearch,
    SearchResult,
)

logger = logging.getLogger("search-tests")


class OmaGroupSearchTests(TestWithDbInstance):
    query = 17
    exp_group = 17

    def test_search_grp_nr(self):
        s = OmaGroupSearch(self.db, self.query)
        self.assertEqual([self.exp_group], [g.group_nbr for g in s.search_groups()])

    def test_search_grp_proteins(self):
        s = OmaGroupSearch(self.db, self.query)
        search_res = s.search_entries()
        for p in search_res:
            self.assertEqual(self.exp_group, p.oma_group)


class OmaGroupSearchWithFingerprintTests(OmaGroupSearchTests):
    query = "HAISGRE"


class GOSearchTest(TestWithDbInstance):
    def test_go_term_search(self):
        query = "GO:0000501"
        s = GOSearch(self.db, query)
        res = [e.omaid for e in s.search_entries()]
        self.assertIn("YEAST00011", res)

    def test_invalid_term_does_list_nothing(self):
        s = GOSearch(self.db, "GO:000000")
        self.assertEqual(0, len(s.search_entries()))


class ECSearchTest(TestWithDbInstance):
    def test_search_full_ec_term(self):
        query = "5.21.5.7"
        with patch.object(self.db, "entrynrs_with_ec_annotation") as mocked:
            mocked.return_value = [12, 4364]  # testdb has no ec annotations
            s = ECSearch(self.db, query)
            self.assertEqual([12, 4364], [p.entry_nr for p in s.search_entries()])


class TaxSearchTest(TestWithDbInstance):
    def test_existing_tax_of_internal_node(self):
        for query in ("Ascomycota", 4890, "Ascomicotta"):
            with self.subTest("existing internal node", query=query):
                s = TaxSearch(self.db, query)
                self.assertIn(
                    4890, [z.ncbi_taxon_id for z in s.search_ancestral_genomes()]
                )
                self.assertIn(559292, [z.ncbi_taxon_id for z in s.search_species()])

    def test_existing_extant_node(self):
        for query in ("SCHPO", 284812, "Schizosaccharomyces pombe"):
            with self.subTest("existing extant node", query=query):
                s = TaxSearch(self.db, query)
                self.assertIn(284812, [z.ncbi_taxon_id for z in s.search_species()])

    def test_as_entries_contains_right_ranges(self):
        query = "Ascomycota"
        s = TaxSearch(self.db, query)
        expected_entries = set([])
        for g in s.search_species():
            a, b = self.db.id_mapper["OMA"].genome_range(g.uniprot_species_code)
            expected_entries.update({nr for nr in range(a, b + 1)})
        self.assertSetEqual(expected_entries, s.search_entries())


class SeqSearchTest(TestWithDbInstance):
    def test_search_example(self):
        query = "ADRIAN"
        s = SequenceSearch(self.db, query)
        res = s.search_entries()
        self.assertGreaterEqual(len(res), 1)
        for p in res:
            self.assertIn(query, p.sequence)

    def test_uniprot_formated_seq(self):
        query = """        10         20         30         40         50
MVKETKFYDI LGVPVTATDV EIKKAYRKCA LKYHPDKNPS EEAAEKFKEA
        60         70         80         90        100
SAAYEILSDP EKRDIYDQFG EDGLSGAGGA GGFPGGGFGF GDDIFSQFFG
       110        120        130        140        150
AGGAQRPRGP QRGKDIKHEI SASLEELYKG RTAKLALNKQ ILCKECEGRG
"""
        s = SequenceSearch(self.db, query)
        self.assertEqual(len(s.seq), 150)
        res = s.search_entries()
        self.assertGreaterEqual(len(res), 1)
        for p in res:
            self.assertIn(s.seq.decode(), p.sequence)

    def test_with_existing_fingerprint_as_sequence(self):
        query = "HAISGRE"
        s = SequenceSearch(self.db, query)
        self.assertIn(17, [g.group_nbr for g in s.search_groups()])


class XRefSearchTest(TestWithDbInstance):
    def test_existing_example(self):
        for query in ("YAS", "K", "YEAST12", "14561"):
            with self.subTest(query=query):
                s = XRefSearch(self.db, query)
                hits = 0
                for p in s.search_entries():
                    for source, val in p.xref_data.items():
                        for el in val:
                            self.assertIn(query.lower(), el.lower())
                            hits += 1
                self.assertGreater(hits, 0)

    def test_total_and_limit(self):
        for query in ("K", "Q1"):
            with self.subTest(query=query):
                s_ref = XRefSearch(self.db, query)
                self.assertGreater(s_ref.estimated_occurrences, 100)
                t0 = time.time()
                res_ref = s_ref.search_entries()
                t1 = time.time()
                s_limit = XRefSearch(self.db, query, max_matches=20)
                res_limit = s_limit.search_entries()
                t2 = time.time()
                self.assertEqual(20, len(res_limit))
                self.assertGreater(len(res_ref), len(res_limit))
                self.assertGreater(t1 - t0, t2 - t1, "limited search took longer")


class HogIDSearchTest(TestWithDbInstance):
    def test_existing_hog_with_level(self):
        for query in ("HOG:0000002_4890", "HOG:0000165.1a_4890"):
            with self.subTest(query=query):
                s = HogIDSearch(self.db, query)
                self.assertIn(
                    query.split("_")[0], [h.hog_id for h in s.search_groups()]
                )
                self.assertFalse(s.outdated_query_hog)

    def test_inexact_hogid_with_level(self):
        for query, (exp_hog_id, exp_level, is_root) in zip(
            (
                "HOG:0000002.5a.6c_4890",
                "HOG:0000165.1a_4751",
                "HOG:0000165_4890",
                "HOG:0000165_4751",
            ),
            (
                ("HOG:0000002", "Ascomycota", False),
                ("HOG:0000165.1a", "Ascomycota", True),
                ("HOG:0000165", "Eukaryota", True),
                ("HOG:0000165", "Fungi", False),
            ),
        ):
            with self.subTest(query=query):
                s = HogIDSearch(self.db, query)
                res = s.search_groups()
                self.assertEqual(1, len(res))
                res = res[0]
                self.assertFalse(s.outdated_query_hog)
                self.assertEqual(exp_hog_id, res.hog_id)
                self.assertEqual(exp_level, res.level)
                self.assertEqual(is_root, res.is_root)

    def test_entries_from_hogid(self):
        query = "HOG:0000165.1a"
        s = HogIDSearch(self.db, query)
        self.assertIsNone(s.search_entries())
        tax = TaxSearch(self.db, "Saccharomycetes")
        s.set_taxon_filter(tax.search_ancestral_genomes()[0])
        for e in s.search_entries():
            self.assertIsInstance(e, ProteinEntry)

    def test_bogous_hogid(self):
        query = "HOG:1111111.6a"
        s = HogIDSearch(self.db, query)
        self.assertEqual([], s.search_groups())


class CombineTest(TestWithDbInstance):
    def test_combine_tax_limit(self):
        s1 = TaxSearch(self.db, "Saccharomycetes")
        s2 = SequenceSearch(self.db, "HAISGRE")
        res = SearchResult()
        res &= s1
        self.assertEqual(len(s1.search_entries()), len(res.entries_set))
        self.assertGreater(len(s2.search_entries()), 2)
        res &= s2
        self.assertEqual(len(s2.search_entries()), len(res.entries_set))
        self.assertEqual(2, len(res.species))
        self.assertEqual(1, len(res.ancestral_genomes))
        self.assertEqual(1, len(res.groups))

    def test_combine_tax_and_hog(self):
        s1 = TaxSearch(self.db, "Saccharomycetes")
        s2 = HogIDSearch(self.db, "HOG:0000165")
        res = SearchResult() & s1 & s2
        self.assertTrue(
            any(
                x.hog_id.startswith("HOG:0000165.1a") and x.level == "Saccharomycetes"
                for x in res.groups.values()
                if isinstance(x, HOG)
            ),
            str(res.groups),
        )

    def test_combine_tax_and_xref(self):
        s1 = TaxSearch(self.db, "Saccharomycetes")
        s2 = XRefSearch(self.db, "K")
        unfiltered_xref = s2.search_entries()
        res = SearchResult() & s1 & s2
        self.assertIsNotNone(s2.entry_filter)
        tax_range = s1.search_entries()
        self.assertTrue(all(p.entry_nr in tax_range for p in res.entries.values()))
        self.assertGreater(len(unfiltered_xref), len(res.entries))
