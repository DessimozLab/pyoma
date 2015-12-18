import unittest
import re
import collections
import numpy
from pyoma.browser.db import *


class TestHelperFunctions(unittest.TestCase):
    def test_counter(self):
        self.assertEqual(0, count_elements([]))
        self.assertEqual(3, count_elements('abc'))
        recarray = numpy.zeros(2, dtype=[('A','i4'),('B','f8')])
        self.assertEqual(2, count_elements(recarray))


class DatabaseTests(unittest.TestCase):
    def setUpClass(cls):

        cls.db = Database()

