from __future__ import absolute_import, unicode_literals, print_function, division
import unittest
import tempfile
import os
import shutil
import lxml.etree

try:
    import unittest.mock as mock
except ImportError:
    import mock as mock
from pyoma.browser.OrthoXMLSplitter import OrthoXMLSplitter


TEST_ORTHXML = """<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="Test" originVersion="0.2">
<species name="HUMAN" NCBITaxId="9601">
  <database name="HUMANfake" version="0.1">
   <genes>
    <gene id="1" protId="HUMAN1" geneId="HUMANg1" />
    <gene id="2" protId="HUMAN2" geneId="HUMANg2" />
    <gene id="3" protId="HUMAN3" geneId="HUMANg3" />
   </genes>
  </database>
 </species>
 <species name="PANTR" NCBITaxId="9483">
  <database name="PANTRfake" version="0.1">
   <genes>
    <gene id="11" protId="PANTR1" geneId="PANTRg1" />
    <gene id="12" protId="PANTR2" geneId="PANTRg2" />
    <gene id="13" protId="PANTR3" geneId="PANTRg3" />
   </genes>
  </database>
 </species>
 <groups>
  <orthologGroup id="1">
      <property name="TaxRange" value="Primates"/>
      <geneRef id="1" />
      <geneRef id="11" />
  </orthologGroup>
  <orthologGroup id="2">
      <property name="TaxRange" value="Primates"/>
      <geneRef id="2" />
      <geneRef id="12" />
  </orthologGroup>
  <orthologGroup id="3">
      <property name="TaxRange" value="Primates"/>
      <geneRef id="3" />
      <geneRef id="13" />
  </orthologGroup>
 </groups>
</orthoXML>
"""


class OrthoXMLSplitterTester(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.infile = os.path.join(cls.tmpdir, "input.orthoxml")
        with open(cls.infile, "w") as fh:
            fh.write(TEST_ORTHXML)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)

    def setUp(self):
        self.outdir = os.path.join(self.tmpdir, "splits")
        self.splitter = OrthoXMLSplitter(self.infile, self.outdir)
        self.splitter.create_new_orthoxml = mock.MagicMock()

    def test_split_into_single_files(self):
        self.splitter()
        self.assertEqual(3, self.splitter.create_new_orthoxml.call_count)

    def test_extract_subset(self):
        self.splitter(hogs_to_extract=[1, 2])
        self.assertEqual(2, self.splitter.create_new_orthoxml.call_count)

    def test_extract_subset_into_single_file(self):
        self.splitter(
            hogs_to_extract=[1, 2], single_hog_files=True, basename="single.orthoxml"
        )
        self.assertEqual(1, self.splitter.create_new_orthoxml.call_count)
        args, kwargs = self.splitter.create_new_orthoxml.call_args
        self.assertEqual(args[0], os.path.join(self.outdir, "single.orthoxml"))

    def test_raises_exception_on_single_hog_files_without_filename(self):
        with self.assertRaises(ValueError):
            self.splitter(hogs_to_extract=[1, 2], single_hog_files=True)


class OrthoXMLSplitterResultTester(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.infile = os.path.join(cls.tmpdir, "input.orthoxml")
        with open(cls.infile, "w") as fh:
            fh.write(TEST_ORTHXML)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)

    def setUp(self):
        self.outdir = os.path.join(self.tmpdir, "splits")
        self.splitter = OrthoXMLSplitter(self.infile, self.outdir)()

    def load_data_of_file(self, fn):
        xml = lxml.etree.parse(fn)
        genes = [
            g.get("id")
            for g in xml.getroot().findall(".//{http://orthoXML.org/2011/}gene")
        ]
        return genes

    def test_properly_split_in_hogs(self):
        for nr in range(1, 4):
            fn = os.path.join(self.outdir, "HOG{:07d}.orthoxml".format(nr))
            exp = [str(nr), str(nr + 10)]
            self.assertEqual(
                sorted(self.load_data_of_file(fn)), sorted(exp), "{} failed.".format(fn)
            )
