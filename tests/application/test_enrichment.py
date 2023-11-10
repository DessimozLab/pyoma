import os
import unittest
import random
from ..browser.test_db import TestWithDbInstance
import pyoma.application.enrichment


class GOEnrichmentTest(TestWithDbInstance):
    def test_yeast_mapk(self):
        res = pyoma.application.enrichment.extant_species_go_enrichment(self.db, [2834, 678, 1003])
        self.assertIn("GO:0000168", res["GO_ID"].tolist())

    def test_random_set(self):
        sample = random.sample(range(*self.db.id_mapper["OMA"].genome_range(1)), k=15)
        res = pyoma.application.enrichment.extant_species_go_enrichment(self.db, sample)
        self.assertEqual(0, len(res[res["p_bonferroni"] < 0.05]))
