from __future__ import unicode_literals, division, absolute_import

import sys
import time
import unittest
from builtins import str

import numpy
from future.utils import with_metaclass

import pyoma.browser.exceptions
from pyoma.browser import models, db, tablefmt
from .test_db import find_path_to_test_db


class TestDbBase(unittest.TestCase):
    db = None

    @classmethod
    def setUpClass(cls):
        fn = find_path_to_test_db("TestDb.h5")
        cls.db = db.Database(fn)

    @classmethod
    def tearDownClass(cls):
        cls.db.close()


class ProteinEntryTests(TestDbBase):
    def test_init_from_enr(self):
        protein_entry = models.ProteinEntry.from_entry_nr(self.db, 12)
        self.assertEqual(protein_entry.entry_nr, 12)
        try:
            protein_entry.canonicalid
        except AttributeError:
            self.assertTrue(False, "entry not properly loaded")

    def test_lazy_eval_of_properties(self):
        protein_entry = models.ProteinEntry.from_entry_nr(self.db, 12)
        self.assertNotIn("sequence", protein_entry.__dict__)
        seq = protein_entry.sequence
        self.assertIn("sequence", protein_entry.__dict__)

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


class HOGModelTest(TestDbBase):
    def test_instant_with_fam(self):
        hog = models.HOG(self.db, 2)
        self.assertEqual(self.db.format_hogid(2), hog.hog_id)
        self.assertTrue(hog.is_root)

    def test_instance_from_hogid(self):
        hog = models.HOG(self.db, self.db.format_hogid(1))
        self.assertEqual(1, hog.fam)
        self.assertTrue(hog.is_root)

    def test_invalid_level(self):
        with self.assertRaises(ValueError):
            hog = models.HOG(self.db, 2, "Metazoa")
            hog.hog_id

    def test_with_valid_level(self):
        hog = models.HOG(self.db, 2, "Fungi")
        self.assertEqual("Fungi", hog.level)
        self.assertFalse(hog.is_root)

    def test_members_and_nr_members(self):
        hog = models.HOG(self.db, "HOG:0000002.1a")
        self.assertEqual(2, hog.fam)
        self.assertEqual(hog.nr_member_genes, len(hog.members))

    def test_keyword_of_hog(self):
        hog = models.HOG(self.db, 2)
        self.assertEqual("", hog.keyword)


class OmaGroupModelTest(TestDbBase):
    def test_instant_with_grpnr(self):
        grp = models.OmaGroup(self.db, 1384)
        self.assertEqual(grp.group_nbr, 1384)

    def test_raises_invalid_id_for_invalid_group_nrs(self):
        for invalid_grp_nr in (-1, 0, 5421):
            with self.assertRaises(pyoma.browser.exceptions.InvalidId):
                grp = models.OmaGroup(self.db, invalid_grp_nr)
                grp.group_nbr

    def test_instant_with_fingerprint(self):
        grp = models.OmaGroup(self.db, "SDNEIRR")
        self.assertEqual("SDNEIRR", grp.fingerprint)

    def test_grp_size(self):
        for grpnr in (521, 4125, 532, 12):
            grp = models.OmaGroup(self.db, grpnr)
            self.assertEqual(len(grp), len(grp.members))

    def test_keyword(self):
        grp = models.OmaGroup(self.db, 521)
        self.assertEqual("methylenetetrahydrofolate reductase", grp.keyword)


class GenomeModelTest(TestDbBase):
    def test_count_genes(self):
        g = self.db.id_mapper["OMA"].identify_genome("YEAST")
        genome = models.Genome(self.db, g)
        self.assertGreaterEqual(genome.nr_genes, genome.nr_entries)

    def test_nr_genes_adds_up(self):
        genomes = self.db.id_mapper["OMA"].genome_table
        tot_genes = sum(models.Genome(self.db, g).nr_genes for g in genomes)
        self.assertEqual(tot_genes, self.db.count_main_isoforms())

    def test_chromosome_length(self):
        g = self.db.id_mapper["OMA"].identify_genome("YEAST")
        genome = models.Genome(self.db, g)
        chr_len = genome.approx_chromosome_length("I")
        self.assertGreaterEqual(chr_len, 220000)
        self.assertLessEqual(chr_len, 230218)  # actual len according to Ensembl

    def test_instantiate_from_species_code(self):
        genome = models.Genome(self.db, "YEAST")
        self.assertEqual(genome.uniprot_species_code, "YEAST")

    def test_raises_for_invalid_species(self):
        with self.assertRaises(pyoma.browser.exceptions.UnknownSpecies):
            genome = models.Genome(self.db, "HUMAN")
            genome.uniprot_species_code

    def test_instantitate_from_taxid(self):
        genome = models.Genome(self.db, 559292)
        self.assertEqual(genome.uniprot_species_code, "YEAST")


class AncestralGenomeModelTests(TestDbBase):
    def test_ncbi_taxon_id_agrees(self):
        query = 4890
        ag = models.AncestralGenome(self.db, query)
        self.assertEqual(query, ag.ncbi_taxon_id)

    def test_retrieve_extant_genomes(self):
        query = 4890
        ag = models.AncestralGenome(self.db, query)
        self.assertIn("YEAST", (g.uniprot_species_code for g in ag.extant_genomes))

    def test_lineage(self):
        query = 4890
        ag = models.AncestralGenome(self.db, query)
        self.assertEqual(ag.lineage[0], ag.sciname)
        self.assertEqual("Eukaryota", ag.kingdom)

    def test_raises_UnknownSpecies_for_invalid_taxid(self):
        with self.assertRaises(pyoma.browser.exceptions.UnknownSpecies):
            tax66 = models.AncestralGenome(self.db, 66).ncbi_taxon_id

    def test_with_luca(self):
        query = "LUCA"
        ag = models.AncestralGenome(self.db, query)
        self.assertEqual(0, ag.ncbi_taxon_id)

    def test_extant_genomes_of_luca(self):
        self.assertEqual(4, len(models.AncestralGenome(self.db, 0).extant_genomes))

    def test_extant_genome_of_extant_genome_taxid(self):
        query = 559292
        ag = models.AncestralGenome(self.db, query)
        self.assertIn(query, [g.ncbi_taxon_id for g in ag.extant_genomes])

    def test_nr_genes(self):
        query_level = "Eukaryota"
        ag = models.AncestralGenome(self.db, query_level)
        self.assertLess(100, ag.nr_genes)


class PairwiseRelationModelTest(TestDbBase):
    def test_relations(self):
        tab_data = self.db.get_vpairs_between_species_pair("YEAST", "PLAF7")
        rels = [models.PairwiseRelation(self.db, d) for d in tab_data]
        self.assertEqual(tab_data[0]["EntryNr1"], rels[0].entry_1.entry_nr)
        self.assertEqual(tab_data[0]["Distance"], rels[0].distance)

    def test_faster_when_passing_entries(self):
        tab_data = self.db.get_vpairs_between_species_pair("YEAST", "PLAF7")
        entries = {e["EntryNr"]: models.ProteinEntry(self.db, e) for e in self.db.main_isoforms("YEAST")}

        t0 = time.time()
        rels = [models.PairwiseRelation(self.db, d) for d in tab_data]
        entry1_nrs_1 = [r.entry_1.entry_nr for r in rels]
        t1 = time.time()
        rels2 = [models.PairwiseRelation(self.db, d, entry1=entries[d["EntryNr1"]]) for d in tab_data]
        entry1_nrs_2 = [r.entry_1.entry_nr for r in rels2]
        t2 = time.time()
        self.assertEqual(entry1_nrs_1, entry1_nrs_2)
        print(t2 - t1, "vs", t1 - t0, "seconds elapsed")
        self.assertLessEqual(t2 - t1, t1 - t0)


class ExonStructureTest(unittest.TestCase):
    def get_exons(self):
        loc_dtype = tablefmt.tables.dtype_from_descr(tablefmt.LocusTable)
        return numpy.array([(1, 500, 510, 1), (1, 600, 610, 1)], dtype=loc_dtype)

    def test_exon_struct_len(self):
        ex = models.ExonStructure(None, self.get_exons())
        self.assertEqual(2, len(ex))

    def test_str_repr(self):
        ex = models.ExonStructure(None, self.get_exons())
        self.assertEqual("join(500..510, 600..610)", str(ex))

    def test_json_repr(self):
        ex = models.ExonStructure(None, self.get_exons())
        self.assertEqual(
            [
                {"start": 500, "end": 510, "strand": "+"},
                {"start": 600, "end": 610, "strand": "+"},
            ],
            ex.as_list_of_dict(),
        )

    def test_str_repr_if_reverse_complement(self):
        ex_dat = self.get_exons()
        ex_dat["Strand"] = -1
        ex = models.ExonStructure(None, ex_dat)
        self.assertEqual("join(complement(600..610), complement(500..510))", str(ex))


class GeneOntologyAnnotationTests(unittest.TestCase):
    db = None

    @classmethod
    def setUpClass(cls):
        fn = find_path_to_test_db("TestDb.h5")
        cls.db = db.Database(fn)

    @classmethod
    def tearDownClass(cls):
        cls.db.get_hdf5_handle().close()

    def test_from_annotation_table(self):
        annos = self.db.get_gene_ontology_annotations(12)
        goa = models.GeneOntologyAnnotation(self.db, annos[0])
        self.assertEqual(12, goa.entry_nr)
        self.assertIn(
            goa.aspect,
            ["molecular_function", "biological_process", "cellular_component"],
        )
        self.assertEqual("YEAST00012", goa.object_id)


class SingletonTests(unittest.TestCase):
    def test_singleton(self):
        class Foo(with_metaclass(models.Singleton, object)):
            def __init__(self):
                pass

        a = Foo()
        b = Foo()
        self.assertEqual(a, b)
