import os
import pickle
import random
import subprocess
import sys
import tempfile
import unittest

import ete4
import numpy
import tables

from pyoma.browser.hogprofile.species_tree import (
    clade_ranges,
    ete_tree_to_rows,
    load_legacy_pickled_tree,
    load_species_tree,
    store_species_tree,
    taxid_from_node_name,
    topmost_nodes_matching,
    tree_to_rows,
)
from pyoma.browser.hogprofile.tree_helper import leaf_index
from .test_db import TestWithDbInstance, find_path_to_test_db


def nested_tree_to_rows(tree):
    """tree given as nested (name, [children]) tuples"""
    return tree_to_rows(tree, children_of=lambda n: n[1], name_of=lambda n: n[0])


def random_newick(nr_leaves, seed=1):
    rng = random.Random(seed)
    nodes = [str(i) for i in range(1, nr_leaves + 1)]
    internal = nr_leaves + 1
    while len(nodes) > 1:
        k = min(len(nodes), rng.randint(1, 4))
        idx = sorted(rng.sample(range(len(nodes)), k), reverse=True)
        children = [nodes.pop(i) for i in idx]
        name = "{}Rep".format(internal) if rng.random() < 0.1 else str(internal)
        nodes.append("({}){}".format(",".join(children), name))
        internal += 1
    return nodes[0] + ";"


class TaxidFromNodeNameTest(unittest.TestCase):
    def test_plain_rep_and_wrapped_uint32(self):
        self.assertEqual(2759, taxid_from_node_name("2759"))
        self.assertEqual(2759, taxid_from_node_name("2759Rep"))
        self.assertEqual(-96, taxid_from_node_name(str(2**32 - 96)))
        self.assertEqual(-5, taxid_from_node_name(-5))


class TreeToRowsTest(unittest.TestCase):
    def setUp(self):
        #        r(1)
        #   x(2)        z(3)
        #  a(4) b(5)  c(6)  y(7Rep)
        #                  d(8) e(9)
        tree = ("1", [("2", [("4", []), ("5", [])]), ("3", [("6", []), ("7Rep", [("8", []), ("9", [])])])])
        self.rows = nested_tree_to_rows(tree)

    def test_preorder_rows(self):
        self.assertEqual([1, 2, 4, 5, 3, 6, 7, 8, 9], self.rows["TaxonId"].tolist())
        self.assertEqual([-1, 0, 1, 1, 0, 4, 4, 6, 6], self.rows["ParentIdx"].tolist())
        self.assertEqual([-1, -1, 0, 1, -1, 2, -1, 3, 4], self.rows["LeafIdx"].tolist())

    def test_clade_ranges(self):
        start, stop, size = clade_ranges(self.rows)
        self.assertEqual([0, 0, 0, 1, 2, 2, 3, 3, 4], start.tolist())
        self.assertEqual([5, 2, 1, 2, 5, 3, 5, 4, 5], stop.tolist())
        self.assertEqual([5, 2, 1, 1, 3, 1, 2, 1, 1], size.tolist())

    def test_topmost_nodes_matching_stops_descending(self):
        mask = numpy.isin(self.rows["TaxonId"], [2, 6, 7, 8])
        self.assertEqual([1, 5, 6], topmost_nodes_matching(self.rows, mask).tolist())
        mask = numpy.isin(self.rows["TaxonId"], [1, 2])
        self.assertEqual([0], topmost_nodes_matching(self.rows, mask).tolist())


class CompareWithEteTreeTest(unittest.TestCase):
    """compares against the tree based implementation of the previous Profiler"""

    def setUp(self):
        self.ete = ete4.PhyloTree(random_newick(400), parser=1)
        self.rows = ete_tree_to_rows(self.ete)

    def reference_annotation(self):
        """logic of the previous tree based Profiler"""
        leaf_idx = leaf_index(self.ete)

        def idx_map(node):
            if node.is_leaf:
                idx = leaf_idx[node.name]
                node.add_props(range=(idx, idx + 1), size=1)
            else:
                min_range, max_range, size = 10000, 0, 0
                for child in node.children:
                    idx_map(child)
                    min_range = min(min_range, child.props["range"][0])
                    max_range = max(max_range, child.props["range"][1])
                    size += child.props["size"]
                node.add_props(range=(min_range, max_range), size=size)

        idx_map(self.ete)
        return leaf_idx

    def test_leaf_idx_matches_profile_columns_of_builder(self):
        leaf_idx = self.reference_annotation()
        leaves = self.rows[self.rows["LeafIdx"] >= 0]
        self.assertEqual({int(k): v for k, v in leaf_idx.items()}, dict(zip(leaves["TaxonId"], leaves["LeafIdx"])))

    def test_clade_ranges_as_tree(self):
        self.reference_annotation()
        start, stop, size = clade_ranges(self.rows)
        nodes = list(self.ete.traverse("preorder"))
        self.assertEqual([n.props["range"] for n in nodes], list(zip(start.tolist(), stop.tolist())))
        self.assertEqual([n.props["size"] for n in nodes], size.tolist())

    def test_key_level_selection_as_tree(self):
        self.reference_annotation()
        _, _, sizes = clade_ranges(self.rows)
        for size in (10, 50, 200):
            expected = [
                taxid_from_node_name(n.name)
                for n in self.ete.leaves(is_leaf_fn=lambda n: size >= n.props["size"] > 0.1 * size)
            ]
            self.assertGreater(len(expected), 1)
            selected = topmost_nodes_matching(self.rows, (sizes <= size) & (sizes > 0.1 * size))
            self.assertEqual(expected, self.rows["TaxonId"][selected].tolist())


class LegacyPickledTreeTest(unittest.TestCase):
    def setUp(self):
        self.db_path = find_path_to_test_db("TestDb.h5")

    def test_load_without_ete3(self):
        # simulate an environment without ete3 installed
        code = (
            "import sys; sys.modules['ete3'] = None\n"
            "import tables\n"
            "from pyoma.browser.hogprofile.species_tree import load_species_tree\n"
            "with tables.open_file(sys.argv[1]) as h5:\n"
            "    rows = load_species_tree(h5.get_node('/HOGProfile/ALL'))\n"
            "print(','.join(map(str, rows['TaxonId'])))\n"
        )
        res = subprocess.run([sys.executable, "-c", code, self.db_path], capture_output=True, text=True)
        self.assertEqual(0, res.returncode, res.stderr)
        taxids = res.stdout.strip().split(",")
        self.assertEqual("2759", taxids[0])
        self.assertEqual(17, len(taxids))

    def test_rows_of_legacy_tree(self):
        # expected rows were verified against the tree unpickled with ete3
        with tables.open_file(self.db_path) as h5:
            rows = load_species_tree(h5.get_node("/HOGProfile/ALL"))
        expected = [
            (2759, -1, -1), (33154, 0, -1), (4751, 1, -1), (451864, 2, -1), (4890, 3, -1), (451866, 4, -1),
            (284812, 5, 0), (716545, 4, -1), (147537, 7, -1), (4891, 8, -1), (4892, 9, -1), (4893, 10, -1),
            (33170, 11, -1), (284811, 12, 1), (559292, 11, 2), (33630, 0, -1), (36329, 15, 3),
        ]  # fmt: skip
        self.assertEqual(expected, rows.tolist())

    def test_rejects_unexpected_globals(self):
        with self.assertRaises(pickle.UnpicklingError):
            load_legacy_pickled_tree(pickle.dumps(tempfile.TemporaryFile))


class StoredTableTest(unittest.TestCase):
    def test_roundtrip_gives_same_rows_as_legacy_format(self):
        with tables.open_file(find_path_to_test_db("TestDb.h5")) as h5:
            legacy_rows = load_species_tree(h5.get_node("/HOGProfile/ALL"))
        with tempfile.TemporaryDirectory() as tmpdir:
            fn = os.path.join(tmpdir, "profile.h5")
            with tables.open_file(fn, "w") as h5:
                group = h5.create_group("/HOGProfile", "ALL", createparents=True)
                store_species_tree(h5, group, legacy_rows)
            with tables.open_file(fn) as h5:
                rows = load_species_tree(h5.get_node("/HOGProfile/ALL"))
        numpy.testing.assert_array_equal(legacy_rows, rows)


class ProfilerTest(TestWithDbInstance):
    def test_species_names_and_key_levels(self):
        from pyoma.browser.hogprofile import Profiler

        p = Profiler(self.db)
        self.assertEqual(
            [
                "Schizosaccharomyces pombe (strain 972 / ATCC 24843)",
                "Ashbya gossypii (strain ATCC 10895 / CBS 109.51 / FGSC 9923 / NRRL Y-1056)",
                "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)",
                "Plasmodium falciparum (isolate 3D7)",
            ],
            p.species_names,
        )
        self.assertEqual(
            {50: {}, 200: {}, "root": {"Opisthokonta": (0, 3), "Alveolata": (3, 4)}},
            p.tax_range_index,
        )

    def test_similar_hog_profile_query(self):
        res = self.db.get_families_with_similar_hog_profile(1)
        self.assertIsNotNone(res)
        self.assertEqual(4, len(res.query_profile))
        self.assertEqual(4, len(res.species_names))
