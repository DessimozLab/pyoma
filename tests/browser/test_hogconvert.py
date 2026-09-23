from __future__ import division, print_function, unicode_literals

import copy
import io
import os
import shutil
import tempfile
import unittest

import numpy
import tables
import lxml.etree as etree
from ete3 import TreeNode

from pyoma.browser import tablefmt
from pyoma.browser.build.hogconvert import (
    Annotator,
    Gene,
    GeneLookupHelper,
    HOGtoHDF5,
    HogObserver,
    OrthoXMLGeneIdParserGeneralProtID,
    OrthoXMLGeneIdParserWithOmaProtId,
    PerFamilyHOGObserver,
    Species,
    TaxonomyLookupHelper,
    parse_orthoxml,
    strip_namespace,
)


class StripNamespaceTest(unittest.TestCase):
    def test_strips_namespace_recursively_and_keeps_attribs_and_text(self):
        xml = b'<geneRef xmlns="http://orthoXML.org/2011/" id="5" LOFT="x">' b'<child a="1">hello</child></geneRef>'
        el = etree.fromstring(xml)
        self.assertEqual("{http://orthoXML.org/2011/}geneRef", el.tag)

        stripped = strip_namespace(el)
        self.assertEqual("geneRef", stripped.tag)
        self.assertEqual({"id": "5", "LOFT": "x"}, dict(stripped.attrib))
        self.assertEqual(1, len(stripped))
        child = stripped[0]
        self.assertEqual("child", child.tag)
        self.assertEqual({"a": "1"}, dict(child.attrib))
        self.assertEqual("hello", child.text)


class FakeParser:
    """minimal stand-in for AbstractOrthoXMLParser, only implementing the
    plumbing Annotator.__init__ needs."""

    def __init__(self, rel_char=""):
        self.rel_char = rel_char
        self._annotators = []

    def register_group_annotator(self, a):
        self._annotators.append(a)


class AnnotatorFormatHogidTest(unittest.TestCase):
    def test_no_rel_char(self):
        a = Annotator(FakeParser(rel_char=""))
        self.assertEqual("HOG:0000001", a.format_hogid(fam=1))

    def test_with_rel_char(self):
        a = Annotator(FakeParser(rel_char="A"))
        self.assertEqual("HOG:A0000001", a.format_hogid(fam="1"))

    def test_with_subhog(self):
        a = Annotator(FakeParser(rel_char=""))
        self.assertEqual("HOG:0000001.1a", a.format_hogid(fam=1, subhog="1a"))

    def test_with_taxid(self):
        a = Annotator(FakeParser(rel_char=""))
        self.assertEqual("HOG:0000001_9606", a.format_hogid(fam=1, taxid=9606))

    def test_with_subhog_and_taxid(self):
        a = Annotator(FakeParser(rel_char="A"))
        self.assertEqual("HOG:A0000001.1a_9606", a.format_hogid(fam=1, subhog="1a", taxid=9606))


def _build_taxonomy_fixture():
    """Eukaryota -> {Mammalia -> {HUMAN, MOUSE}, YEAST}"""
    root = TreeNode()
    root.add_features(taxid=1, xml_taxonId=1, name="Eukaryota")
    mammalia = root.add_child(name="Mammalia")
    mammalia.add_features(taxid=2, xml_taxonId=2)
    human = mammalia.add_child(name="HUMAN")
    human.add_features(taxid=9606, xml_taxonId=101)
    mouse = mammalia.add_child(name="MOUSE")
    mouse.add_features(taxid=10090, xml_taxonId=102)
    yeast = root.add_child(name="YEAST")
    yeast.add_features(taxid=4932, xml_taxonId=103)
    for n in root.traverse(strategy="preorder"):
        n.add_feature("size", len(n))
    return root, {"root": root, "mammalia": mammalia, "human": human, "mouse": mouse, "yeast": yeast}


class TaxonomyLookupHelperTest(unittest.TestCase):
    def setUp(self):
        self.tree, self.nodes = _build_taxonomy_fixture()
        self.helper = TaxonomyLookupHelper(self.tree)

    def test_get_node_from_taxid(self):
        self.assertIs(self.nodes["human"], self.helper.get_node_from_taxid(9606))
        self.assertIs(self.nodes["root"], self.helper.get_node_from_taxid("1"))
        with self.assertRaises(KeyError):
            self.helper.get_node_from_taxid(424242)

    def test_get_node_from_taxonId(self):
        self.assertIs(self.nodes["human"], self.helper.get_node_from_taxonId(101))
        with self.assertRaises(KeyError):
            self.helper.get_node_from_taxonId(999999)

    def test_get_mrca_single_leaf_shortcut(self):
        res = self.helper.get_mrca([(9606, "taxid")])
        self.assertIs(self.nodes["human"], res)

    def test_get_mrca_two_species_same_family(self):
        res = self.helper.get_mrca([(9606, "taxid"), (10090, "taxid")])
        self.assertIs(self.nodes["mammalia"], res)

    def test_get_mrca_across_root(self):
        res = self.helper.get_mrca([(9606, "taxid"), (4932, "taxid")])
        self.assertIs(self.nodes["root"], res)

    def test_get_mrca_mixed_source_types(self):
        res = self.helper.get_mrca([(101, "xml"), (10090, "taxid")])
        self.assertIs(self.nodes["mammalia"], res)

    def test_levels_between_skips_adjacent(self):
        self.assertEqual([], list(self.helper.levels_between(self.nodes["human"], self.nodes["mammalia"])))

    def test_levels_between_yields_intermediate_nodes(self):
        res = list(self.helper.levels_between(self.nodes["human"], self.nodes["root"]))
        self.assertEqual([self.nodes["mammalia"]], res)

    def test_as_xml_full_tree(self):
        xml = self.helper.as_xml()
        ids_names = {(int(el.get("id")), el.get("name")) for el in xml.iter("taxon")}
        expected = {(1, "Eukaryota"), (2, "Mammalia"), (9606, "HUMAN"), (10090, "MOUSE"), (4932, "YEAST")}
        self.assertEqual(expected, ids_names)

    def test_as_xml_subtree(self):
        xml = self.helper.as_xml(root_tax_node=self.nodes["mammalia"])
        self.assertEqual("taxonomy", xml.tag)
        self.assertEqual("Mammalia", xml[0].get("name"))
        names = {el.get("name") for el in xml.iter("taxon")}
        self.assertEqual({"Mammalia", "HUMAN", "MOUSE"}, names)


class HogFixtureMixin:
    """Builds a small, real HDF5 fixture (Taxonomy/Genome/Protein tables) and
    a matching orthoxml document with a genuine taxonomic gap so that
    Annotator.insert_missing_levels has real work to do:

        Eukaryota(1) -> Metazoa(3) -> Mammalia(2) -> {HUMAN(9606), MOUSE(10090)}
                      -> YEAST(4932)   [species not part of the HOG]

    root orthologGroup is declared at Eukaryota; the two genes' MRCA is
    Mammalia, two taxonomic ranks below the declared root (Metazoa is
    the level strictly in between and gets a synthetic orthologGroup
    inserted for it; then each geneRef itself gets Mammalia + its own
    species level inserted too, since paralogGroup children get expanded
    all the way to species level).
    """

    XML = b"""<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="Test" originVersion="0.1">
<species name="HUMAN" taxonId="9606">
 <database name="HUMANfake" version="1">
  <genes>
   <gene id="1" protId="HUMAN00001"/>
   <gene id="2" protId="HUMAN00002"/>
  </genes>
 </database>
</species>
<species name="MOUSE" taxonId="10090">
 <database name="MOUSEfake" version="1">
  <genes>
   <gene id="11" protId="MOUSE00001"/>
   <gene id="12" protId="MOUSE00002"/>
  </genes>
 </database>
</species>
<taxonomy>
<taxon id="1" name="Eukaryota">
  <taxon id="3" name="Metazoa">
    <taxon id="2" name="Mammalia">
      <taxon id="9606" name="HUMAN"/>
      <taxon id="10090" name="MOUSE"/>
    </taxon>
  </taxon>
  <taxon id="4932" name="YEAST"/>
</taxon>
</taxonomy>
<groups>
<orthologGroup id="1" taxonId="1">
  <property name="TaxRange" value="Eukaryota"/>
  <paralogGroup>
    <geneRef id="1"/>
    <geneRef id="11"/>
  </paralogGroup>
</orthologGroup>
</groups>
</orthoXML>
"""

    def build_h5_fixture(self, path):
        h5 = tables.open_file(path, mode="w")
        h5.set_node_attr("/", "oma_release_char", "")

        tax = numpy.zeros(3, dtype=tables.dtype_from_descr(tablefmt.TaxonomyTable))
        tax[0] = (1, 0, b"Eukaryota", False, numpy.nan)
        tax[1] = (3, 1, b"Metazoa", False, numpy.nan)
        tax[2] = (2, 3, b"Mammalia", False, numpy.nan)
        h5.create_table("/", "Taxonomy", obj=tax)

        gs = numpy.zeros(3, dtype=tables.dtype_from_descr(tablefmt.GenomeTable))
        gs[0]["UniProtSpeciesCode"] = b"HUMAN"
        gs[0]["NCBITaxonId"] = 9606
        gs[0]["EntryOff"] = 0
        gs[0]["TotEntries"] = 2
        gs[0]["SciName"] = b"Homo sapiens"
        gs[1]["UniProtSpeciesCode"] = b"MOUSE"
        gs[1]["NCBITaxonId"] = 10090
        gs[1]["EntryOff"] = 2
        gs[1]["TotEntries"] = 2
        gs[1]["SciName"] = b"Mus musculus"
        gs[2]["UniProtSpeciesCode"] = b"YEAST"
        gs[2]["NCBITaxonId"] = 4932
        gs[2]["EntryOff"] = 4
        gs[2]["TotEntries"] = 0
        gs[2]["SciName"] = b"Saccharomyces cerevisiae"
        h5.create_table("/", "Genome", obj=gs)

        entries = numpy.zeros(4, dtype=tables.dtype_from_descr(tablefmt.ProteinTable))
        entries["EntryNr"] = [1, 2, 3, 4]
        h5.create_group("/", "Protein")
        h5.create_table("/Protein", "Entries", obj=entries)
        h5.close()

    def parse_fixture(self, path):
        parser = OrthoXMLGeneIdParserWithOmaProtId(path)
        Annotator(parser)
        recorded = {"hog": [], "augmented": []}

        class Recorder(HogObserver):
            def process_hog(self, node):
                # process_hog and process_augmented_hog fire on the *same*
                # mutable lxml Element (augmentation mutates it in place),
                # so snapshot a deep copy to see the pre-augmentation state.
                recorded["hog"].append(copy.deepcopy(node))

            def process_augmented_hog(self, node):
                recorded["augmented"].append(node)

        Recorder(parser)
        parse_orthoxml(io.BytesIO(self.XML), parser)
        return parser, recorded


class HogAnnotationEndToEndTest(HogFixtureMixin, unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.path = os.path.join(self.tmpdir, "fixture.h5")
        self.build_h5_fixture(self.path)
        self.parser, self.recorded = self.parse_fixture(self.path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_one_hog_processed(self):
        self.assertEqual(1, len(self.recorded["hog"]))
        self.assertEqual(1, len(self.recorded["augmented"]))

    def test_pre_augmentation_hog_id_and_generef_rewrite(self):
        hog = self.recorded["hog"][0]
        self.assertEqual("HOG:0000001_1", hog.get("id"))
        self.assertEqual("1", hog.get("taxonId"))
        # geneRef ids rewritten from original xml gene ids (1, 11) to real
        # 1-based EntryNr (HUMAN00001 -> 1, MOUSE00001 -> 3)
        generef_ids = [g.get("id") for g in hog.findall(".//geneRef")]
        self.assertEqual(["1", "3"], generef_ids)
        # no LOFT attribute yet before augmentation
        for g in hog.findall(".//geneRef"):
            self.assertIsNone(g.get("LOFT"))

    def test_pre_augmentation_completeness_score(self):
        hog = self.recorded["hog"][0]
        score = hog.find('./score[@id="CompletenessScore"]')
        self.assertIsNotNone(score)
        # 2 covered species (HUMAN, MOUSE) out of Eukaryota's 3 leaves
        # (HUMAN, MOUSE, YEAST)
        self.assertAlmostEqual(2 / 3, float(score.get("value")), places=3)

    def test_augmented_hog_has_inserted_metazoa_level(self):
        aug = self.recorded["augmented"][0]
        # the direct child of the root should now be the synthetic
        # inserted orthologGroup for Metazoa (the level strictly between
        # the declared root level "Eukaryota" and the genes' real MRCA
        # "Mammalia")
        inserted = aug[-1]
        self.assertEqual("orthologGroup", inserted.tag)
        self.assertEqual("3", inserted.get("taxonId"))
        tax_range = inserted.find('./property[@name="TaxRange"]')
        self.assertEqual("Metazoa", tax_range.get("value"))
        nr_genes = inserted.find('./property[@name="NrMemberGenes"]')
        self.assertEqual("2", nr_genes.get("value"))
        score = inserted.find('./score[@id="CompletenessScore"]')
        self.assertAlmostEqual(1.0, float(score.get("value")), places=3)

    def test_augmented_hog_generef_loft_and_species_levels(self):
        aug = self.recorded["augmented"][0]
        generefs = aug.findall(".//geneRef")
        self.assertEqual(2, len(generefs))
        loft_ids = {g.get("id"): g.get("LOFT") for g in generefs}
        self.assertEqual({"1": "HOG:0000001.1a", "3": "HOG:0000001.1b"}, loft_ids)

        # each geneRef should now be wrapped in an inserted species-level
        # orthologGroup (taxonId 9606 / 10090) whose id encodes the same
        # sub-hog id as the geneRef's LOFT value
        for g in generefs:
            species_level = g.getparent()
            self.assertEqual("orthologGroup", species_level.tag)
            self.assertIn(species_level.get("taxonId"), ("9606", "10090"))
            self.assertEqual(f"{g.get('LOFT')}_{species_level.get('taxonId')}", species_level.get("id"))

    def test_generef_ids_match_real_entrynrs(self):
        entries = self.parser._gstab  # sanity: fixture as expected (HUMAN, MOUSE, YEAST)
        self.assertEqual(3, len(entries))
        aug = self.recorded["augmented"][0]
        generef_ids = sorted(int(g.get("id")) for g in aug.findall(".//geneRef"))
        self.assertEqual([1, 3], generef_ids)


class HOGtoHDF5Test(HogFixtureMixin, unittest.TestCase):
    # Regression tests for two compounding bugs, previously unnoticed since
    # this class had 0% test coverage (see git history for the fix commit):
    #
    # (1) HOGtoHDF5.process_augmented_hog computed the "IsRoot" flag for each
    #     HOG level as `bool(ognode.getparent().tag in ("paralogGroup",
    #     "groups"))`. For the family's TOP-LEVEL orthologGroup, `ognode` IS
    #     that root node -- but `AbstractOrthoXMLParser.process_group` always
    #     runs the node through `strip_namespace()` first, which builds a
    #     brand new DETACHED `etree.Element(...)` tree with no parent at all.
    #     So `ognode.getparent()` was always None for the root orthologGroup's
    #     own TaxRange entry, and `.tag` on None raised AttributeError --
    #     unconditionally, for any real orthoxml input, as soon as processing
    #     reached the root level's own TaxRange. Fixed by treating a None
    #     parent as the true root.
    #
    # (2) That AttributeError used to propagate out of parse_orthoxml() while
    #     still inside the `with HOGtoHDF5(...)` block, so
    #     HOGtoHDF5.__exit__(exc_type=AttributeError, ...) ran its data-write
    #     lines (`numpy.stack(list(self.index.values()))` etc.)
    #     UNCONDITIONALLY regardless of exc_type -- and since no rootHOG's
    #     orthoxml had been stored yet (self.index still empty),
    #     `numpy.stack([])` raised `ValueError: need at least one array to
    #     stack`, masking the original AttributeError, and self.h5 was never
    #     closed (leaking the open file handle). Fixed by only running the
    #     write-out logic when exc_type is None, and always closing self.h5
    #     via `finally`.
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmpdir, ignore_errors=True)
        self.path = os.path.join(self.tmpdir, "fixture.h5")
        self.build_h5_fixture(self.path)

    def test_hogtohdf5_writes_leveltab_and_omahog(self):
        parser = OrthoXMLGeneIdParserWithOmaProtId(self.path)
        Annotator(parser)
        out_path = os.path.join(self.tmpdir, "hogs.h5")
        with HOGtoHDF5(parser, out_path) as writer:
            # mirrors real usage in build/main.py: a PerFamilyHOGObserver
            # feeds each (augmented) HOG's serialized orthoxml back into the
            # HOGtoHDF5 writer via these callbacks.
            PerFamilyHOGObserver(parser, writer.store_orthoxml, writer.store_orthoxml_augmented)
            parse_orthoxml(io.BytesIO(self.XML), parser)
        self.assertEqual(0, writer.h5.isopen, "file must be closed after a successful run")

        with tables.open_file(out_path, mode="r") as h5out:
            leveltab = h5out.get_node("/HogLevel").read()
            by_level = {row["Level"].decode(): row for row in leveltab}
            self.assertEqual({"Eukaryota", "Metazoa", "Mammalia", "Homo sapiens", "Mus musculus"}, set(by_level))
            # root of the family: directly under <groups> (detached parent) -> IsRoot
            self.assertTrue(bool(by_level["Eukaryota"]["IsRoot"]))
            self.assertEqual(b"HOG:0000001", by_level["Eukaryota"]["ID"])
            # inserted intermediate level, nested under another orthologGroup -> not root
            self.assertFalse(bool(by_level["Metazoa"]["IsRoot"]))
            # directly under the paralogGroup (start of a new sub-hog lineage) -> IsRoot
            self.assertTrue(bool(by_level["Mammalia"]["IsRoot"]))
            # species-level leaf, nested under the Mammalia orthologGroup -> not root
            self.assertFalse(bool(by_level["Homo sapiens"]["IsRoot"]))
            self.assertFalse(bool(by_level["Mus musculus"]["IsRoot"]))
            self.assertTrue(all(row["Fam"] == 1 for row in leveltab))
            # two rows carry the per-paralog sub-hog id (one per geneRef lineage)
            sub_ids = {row["ID"] for row in leveltab if row["Level"] in (b"Mammalia",)}
            self.assertEqual({b"HOG:0000001.1a", b"HOG:0000001.1b"}, sub_ids)

            oma_hog = h5out.get_node("/OmaHOG").read()
            # entries 1 (HUMAN) and 3 (MOUSE) got their sub-hog id; entries
            # 2, 4 (not part of this HOG) stay empty
            self.assertEqual(b"HOG:0000001.1a", oma_hog[0])
            self.assertEqual(b"", oma_hog[1])
            self.assertEqual(b"HOG:0000001.1b", oma_hog[2])
            self.assertEqual(b"", oma_hog[3])

            index = h5out.get_node("/OrthoXML/Index").read()
            self.assertEqual(1, len(index))
            self.assertEqual(1, index[0]["Fam"])
            self.assertGreater(index[0]["HogAugmentedBufferLength"], 0)

    def test_hogtohdf5_reraises_original_error_and_closes_file(self):
        # a geneRef referencing a gene id that was never registered makes
        # Annotator.update_hogids fail with AttributeError while resolving
        # the gene (`gene_accessor_func(...)` returns None)
        broken_xml = self.XML.replace(b'<geneRef id="1"/>', b'<geneRef id="999"/>', 1)
        parser = OrthoXMLGeneIdParserWithOmaProtId(self.path)
        Annotator(parser)
        out_path = os.path.join(self.tmpdir, "hogs_broken.h5")
        writer_holder = {}
        with self.assertRaises(AttributeError):
            with HOGtoHDF5(parser, out_path) as writer:
                writer_holder["writer"] = writer
                PerFamilyHOGObserver(parser, writer.store_orthoxml, writer.store_orthoxml_augmented)
                parse_orthoxml(io.BytesIO(broken_xml), parser)
        self.assertEqual(0, writer_holder["writer"].h5.isopen, "file must be closed even when processing fails")


class GeneLookupHelperTest(unittest.TestCase):
    def test_standalone_add_species_and_lookup(self):
        gh = GeneLookupHelper(parser=None)
        sp1 = Species("AAAA", 1, 0, 2, xml_node=etree.Element("species"))
        sp2 = Species("BBBB", 2, 2, 5, xml_node=etree.Element("species"))
        sp3 = Species("CCCC", 3, 5, 8, xml_node=etree.Element("species"))
        genes1 = [Gene(1, 1, "orig1", "AAAA0001", 1), Gene(2, 2, "orig2", "AAAA0002", 1)]
        genes2 = [Gene(11, 3, "orig11", "BBBB0001", 2), Gene(12, 5, "orig12", "BBBB0003", 2)]
        genes3 = [Gene(21, 6, "orig21", "CCCC0001", 3)]
        # add out of offset-order to check _sort_and_freeze sorts them
        gh.add_species(sp3, genes3)
        gh.add_species(sp1, genes1)
        gh.add_species(sp2, genes2)

        self.assertEqual(genes1[0], gh.get_gene_by_original_id(1))
        self.assertEqual(genes2[1], gh.get_gene_by_new_id(5))
        self.assertIsNone(gh.get_gene_by_original_id(999))

        # boundaries: first entry of a species, last entry, and the
        # species with the highest offset
        self.assertEqual("AAAA", gh.search_species_by_entry_nr(1).species_code)
        self.assertEqual("AAAA", gh.search_species_by_entry_nr(2).species_code)
        self.assertEqual("BBBB", gh.search_species_by_entry_nr(3).species_code)
        self.assertEqual("BBBB", gh.search_species_by_entry_nr(5).species_code)
        self.assertEqual("CCCC", gh.search_species_by_entry_nr(6).species_code)


class GeneLookupHelperIterSpeciesTest(HogFixtureMixin, unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmpdir, ignore_errors=True)
        self.path = os.path.join(self.tmpdir, "fixture.h5")
        self.build_h5_fixture(self.path)
        self.parser, self.recorded = self.parse_fixture(self.path)

    def test_iter_species_with_genes_nodes_all(self):
        gh = self.parser.gene_helper
        species_nodes = list(gh.iter_species_with_genes_nodes())
        names = [n.get("name") for n in species_nodes]
        self.assertEqual(["Homo sapiens", "Mus musculus"], names)
        human_gene_ids = [g.get("id") for g in species_nodes[0].findall(".//gene")]
        mouse_gene_ids = [g.get("id") for g in species_nodes[1].findall(".//gene")]
        self.assertEqual(["1", "2"], human_gene_ids)
        self.assertEqual(["3", "4"], mouse_gene_ids)

    def test_iter_species_with_genes_nodes_subset(self):
        gh = self.parser.gene_helper
        species_nodes = list(gh.iter_species_with_genes_nodes(new_generef_ids=[1]))
        self.assertEqual(1, len(species_nodes))
        self.assertEqual("Homo sapiens", species_nodes[0].get("name"))
        gene_ids = [g.get("id") for g in species_nodes[0].findall(".//gene")]
        self.assertEqual(["1"], gene_ids)


class SearchEntryNrsOfIdsWithinSpeciesTest(unittest.TestCase):
    def setUp(self):
        self.p = object.__new__(OrthoXMLGeneIdParserGeneralProtID)
        # 5 protein ids at positions 0..4; sp restricted to [0,3)
        self.p._prot_ids = numpy.array([b"P00003", b"P00001", b"P00002", b"P00001", b"Q99999"], dtype="S10")
        self.p._protkey = self.p._prot_ids.argsort()

    def test_resolves_ids_within_range(self):
        sp = {"EntryOff": 0, "TotEntries": 3}
        ids = numpy.array([b"P00002", b"P00003"], dtype="S10")
        res = self.p._search_entry_nrs_of_ids_within_species(ids, sp)
        numpy.testing.assert_array_equal([3, 1], res)

    def test_id_not_present_raises(self):
        sp = {"EntryOff": 0, "TotEntries": 3}
        with self.assertRaises(ValueError):
            self.p._search_entry_nrs_of_ids_within_species(numpy.array([b"ZZZZZZ"], dtype="S10"), sp)

    def test_id_outside_species_range_raises(self):
        # Q99999 exists (position 4) but sp only covers [0,3)
        sp = {"EntryOff": 0, "TotEntries": 3}
        with self.assertRaises(ValueError):
            self.p._search_entry_nrs_of_ids_within_species(numpy.array([b"Q99999"], dtype="S10"), sp)

    def test_duplicate_id_within_range_raises(self):
        # P00001 occurs twice (positions 1 and 3); widen sp to cover both
        sp = {"EntryOff": 0, "TotEntries": 5}
        with self.assertRaises(ValueError):
            self.p._search_entry_nrs_of_ids_within_species(numpy.array([b"P00001"], dtype="S10"), sp)


if __name__ == "__main__":
    unittest.main()
