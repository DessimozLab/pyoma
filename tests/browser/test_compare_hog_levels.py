import unittest
import numpy
from pyoma.browser.hoghelper import compare_levels


class CompareHogLevels(unittest.TestCase):
    def setUp(self) -> None:
        self.parent = numpy.array(
            [
                (1, b"HOG:0000001", b"Parent"),
                (2, b"HOG:0000002", b"Parent"),
                (3, b"HOG:0000003", b"Parent"),
                (4, b"HOG:0000004", b"Parent"),
                (6, b"HOG:0000006.1a", b"Parent"),
                (6, b"HOG:0000006.1b", b"Parent"),
                (7, b"HOG:0000007.2c", b"Parent"),
            ],
            dtype=[("Fam", "i4"), ("ID", "S255"), ("Level", "S255")],
        )
        self.child = numpy.array(
            [
                (1, b"HOG:0000001", b"Child"),
                (2, b"HOG:0000002.1a", b"Child"),
                (4, b"HOG:0000004.1a", b"Child"),
                (4, b"HOG:0000004.1b", b"Child"),
                (5, b"HOG:0000005", b"Child"),
                (6, b"HOG:0000006.1a", b"Child"),
                (7, b"HOG:0000007.2c", b"Child"),
            ],
            dtype=self.parent.dtype,
        )

    def test_full_set_case(self):
        parent = self.parent
        child = self.child
        exp_res = [
            (child[0], "retained"),
            (child[1], "duplicated"),
            (parent[2], "lost"),
            (child[2], "duplicated"),
            (child[3], "duplicated"),
            (child[4], "gained"),
            (child[5], "retained"),
            (parent[5], "lost"),
            (child[6], "retained"),
        ]
        self.assertEqual(exp_res, compare_levels(parent, child))

    def test_fam_1_to_5(self):
        child = self.child[self.child["Fam"] <= 5]
        parent = self.parent[self.parent["Fam"] <= 5]
        exp_res = [
            (child[0], "retained"),
            (child[1], "duplicated"),
            (parent[2], "lost"),
            (child[2], "duplicated"),
            (child[3], "duplicated"),
            (child[4], "gained"),
        ]
        self.assertEqual(exp_res, compare_levels(parent, child))

    def test_fam_2_to_6(self):
        child = self.child[(1 < self.child["Fam"]) & (self.child["Fam"] <= 6)]
        parent = self.parent[(1 < self.child["Fam"]) & (self.parent["Fam"] <= 6)]
        exp_res = [
            (child[0], "duplicated"),
            (parent[1], "lost"),
            (child[1], "duplicated"),
            (child[2], "duplicated"),
            (child[3], "gained"),
            (child[4], "retained"),
            (parent[4], "lost"),
        ]
        self.assertEqual(exp_res, compare_levels(parent, child))


if __name__ == "__main__":
    unittest.main()
