import unittest
from pyoma.browser import hoghelper


class TestAreOrthologsFunction(unittest.TestCase):
    def test_basic_tests(self):
        cases = [
            (("HOG:A0002412", "HOG:A0002412"), True, "HOG:A0002412"),
            (("HOG:A0002412.4a", "HOG:A0002412.4a"), True, "HOG:A0002412.4a"),
            (("HOG:A0002412.4a.6b", "HOG:A0002412.4a"), True, "HOG:A0002412.4a"),
            (("HOG:A0002412.4a", "HOG:A0002412.4a.6b"), True, "HOG:A0002412.4a"),
            (("HOG:A0002412.4a", "HOG:A0002412.5c"), True, "HOG:A0002412"),
            (("HOG:A0001413.4a", "HOG:A0002412"), False, False),
            (("HOG:A0002412.4a", "HOG:A0002412.4b"), False, False),
            (("HOG:A0002412.4a", "HOG:A0002412.4b.6b"), False, False),
            (("HOG:A0002412.14aa", "HOG:A0002412.14ab"), False, False),
            (("HOG:A0002412.15aa", "HOG:A0002412.14aa"), True, "HOG:A0002412"),
            (
                ("HOG:A0002412.15aa.8b", "HOG:A0002412.15aa.3a"),
                True,
                "HOG:A0002412.15aa",
            ),
            (("HOG:A0002412.15a", "HOG:A0002412.14b"), True, "HOG:A0002412"),
            (("HOG:A0002412.1a.4b.2c", "HOG:A0002412.1a.4a.2c"), False, False),
        ]
        for case, exp, grp in cases:
            with self.subTest(case=case, exp=exp, grp=grp):
                self.assertEqual(
                    exp,
                    hoghelper.are_orthologous(*case),
                    f"{case[0]} / {case[1]} wrong orthologous assignment",
                )
                self.assertEqual(
                    grp,
                    hoghelper.are_orthologous(*case, return_group=True),
                    f"{case[0]} / {case[1]} wrong group assignment",
                )
