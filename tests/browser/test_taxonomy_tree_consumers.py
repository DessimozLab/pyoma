"""tests for the code consuming (ete4) trees produced by omataxonomy.Taxonomy.get_topology"""

import unittest

import ete4

from pyoma.browser.build.taxonomy_alignment import collect_ncbi_map_info, create_treenode_to_genome_map
from pyoma.browser.build.xref.filter import load_relevant_taxids


class FakeTaxonomy:
    """mimics the parts of omataxonomy.Taxonomy used by pyoma.

    2759 Eukaryota
    ├── 10 genus G1
    │   ├── 100 species (── 1000 strain)
    │   └── 101 species
    ├── 20 genus G2
    │   └── 200 species
    └── 30 family F3
        └── 31 genus G3
            ├── 300 species
            └── 301 species
    """

    parent = {10: 2759, 100: 10, 1000: 100, 101: 10, 20: 2759, 200: 20, 30: 2759, 31: 30, 300: 31, 301: 31}
    rank = {2759: "superkingdom", 10: "genus", 100: "species", 1000: "strain", 101: "species", 20: "genus"}
    rank.update({200: "species", 30: "family", 31: "genus", 300: "species", 301: "species"})

    def get_lineage(self, taxid):
        lineage = [taxid]
        while lineage[-1] in self.parent:
            lineage.append(self.parent[lineage[-1]])
        return lineage[::-1]

    def get_rank(self, taxids):
        return {t: self.rank[t] for t in taxids if t in self.rank}

    def get_descendant_taxa(self, taxid, intermediate_nodes=False):
        return [t for t in self.parent if t != taxid and taxid in self.get_lineage(t)]

    def get_taxid_translator(self, taxids):
        return {t: f"taxon {t}" for t in taxids if t in self.rank}

    def translate_to_names(self, taxids):
        return [f"taxon {t}" for t in taxids]

    def get_topology(self, taxids, intermediate_nodes=True, annotate=True):
        nodes = {}
        for taxid in taxids:
            for t in self.get_lineage(taxid):
                if t not in nodes:
                    node = ete4.PhyloTree()
                    node.name = str(t)
                    node.add_props(taxid=t, rank=self.rank[t])
                    nodes[t] = node
                    if t in self.parent:
                        nodes[self.parent[t]].add_child(node)
        return nodes[2759]


class LoadRelevantTaxidsTest(unittest.TestCase):
    def test_relevant_taxids_include_subtaxa_and_genus_clades(self):
        res = load_relevant_taxids([100, 101, 1000, 200, 300], FakeTaxonomy())
        expected = {
            100: {100, 1000},
            1000: {100, 1000},
            101: {101},
            200: {200},
            300: {300},
            10: {101, 1000},
            20: {200},
            31: {300},
        }
        self.assertEqual(expected, dict(res))


class CreateTreenodeToGenomeMapTest(unittest.TestCase):
    def test_genome_sets_of_clades(self):
        tax = FakeTaxonomy()
        tree = tax.get_topology([100, 101, 200, 300, 301])
        tax2code = {100: {"AAAAA"}, 101: {"BBBBB"}, 200: {"CCCCC"}, 300: {"DDDDD"}}
        res = create_treenode_to_genome_map(tree, tax2code)
        self.assertEqual({"AAAAA", "BBBBB", "CCCCC", "DDDDD"}, res[2759])
        self.assertEqual({"AAAAA", "BBBBB"}, res[10])
        self.assertEqual({"DDDDD"}, res[31])
        self.assertEqual(set(), res[301])
        self.assertEqual(res[31], tree.search_nodes(taxid=31).__next__().props["spset"])


class CollectNcbiMapInfoTest(unittest.TestCase):
    def test_representatives_of_clades_not_in_oma_taxonomy(self):
        tax = FakeTaxonomy()
        tree = tax.get_topology([100, 101, 200, 300, 301])
        tax2code = {100: {"AAAAA"}, 101: {"BBBBB"}, 200: {"CCCCC"}, 300: {"DDDDD"}, 301: {"EEEEE"}}
        code2tax = {c: {t} for t, codes in tax2code.items() for c in codes}
        pyomatax = {"NCBITaxonId": {2759, 100, 101, 200, 300, 301}}
        cand_map = {2759: [], 10: [], 30: [], 100: []}
        ancestral, extant = collect_ncbi_map_info(cand_map, tax2code, code2tax, tax, pyomatax, tree)
        self.assertEqual([{"ncbi_taxid": 100, "sciname": "taxon 100"}], extant)
        by_id = {a["ncbi_taxid"]: a for a in ancestral}
        self.assertEqual(True, by_id[2759]["valid_oma_tax"])
        self.assertEqual((False, 100), (by_id[10]["valid_oma_tax"], by_id[10]["representative"]))
        self.assertEqual((False, 300), (by_id[30]["valid_oma_tax"], by_id[30]["representative"]))
