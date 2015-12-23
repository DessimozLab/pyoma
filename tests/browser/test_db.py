import unittest
import numpy
import tables
from pyoma.browser.db import *
from pyoma.browser import tablefmt


class TestHelperFunctions(unittest.TestCase):
    def test_counter(self):
        self.assertEqual(0, count_elements([]))
        self.assertEqual(3, count_elements('abc'))
        recarray = numpy.zeros(2, dtype=[('A','i4'),('B','f8')])
        self.assertEqual(2, count_elements(recarray))


class DatabaseTests(unittest.TestCase):
    db = None

    @classmethod
    def setUpClass(cls):
        path = "/pub/projects/cbrg-oma-browser/Test.Jul2014/data/OmaServer.h5"
        cls.db = Database(path)

    @classmethod
    def tearDownClass(cls):
        cls.db.get_hdf5_handle().close()

    def test_get_vpairs_of_entry_with_orthologs(self):
        for entry_nr, exp_vps_cnt in [(12, 3), (1, 0), (4,1)]:
            vps = self.db.get_vpairs(entry_nr)
            self.assertTrue(isinstance(vps, numpy.ndarray))
            self.assertEqual(exp_vps_cnt, len(vps))
            self.assertEqual(['EntryNr1', 'EntryNr2', 'RelType'],
                             sorted(vps.dtype.fields.keys()))

    def test_neighborhood_close_to_boundary(self):
        query, window = 3, 7
        neighbors, idx = self.db.neighbour_genes(query, window)
        self.assertEqual(query - 1, idx)
        self.assertEqual(neighbors['EntryNr'][idx], query)
        expected_entry_nrs = numpy.arange(1, query + window + 1, dtype='i4')
        self.assertTrue(numpy.array_equal(expected_entry_nrs,
                                          neighbors['EntryNr']))

    def test_hog_family(self):
        entry = numpy.zeros(1, dtype=tables.dtype_from_descr(tablefmt.ProteinTable))
        entry['OmaHOG'] = b""
        with self.assertRaises(Singleton):
            self.db.hog_family(entry[0])
        entry['OmaHOG'] = b"HOG:0000523"
        self.assertEqual(523, self.db.hog_family(entry[0]))

    def test_orthoxml(self):
        xml = self.db.get_orthoxml(33).decode()
        # we simply check that orthoxml starts with <?xml and ends with an orthoxml tag
        self.assertTrue(xml.startswith('<?xml '))
        self.assertTrue(xml.endswith('</orthoXML>\n'))

    def test_hog_lex_range(self):
        cases = [(b'HOG:001', (b'HOG:001', b'HOG:002')),
                 (b'HOG:001.1a.2b', (b'HOG:001.1a.2b', b'HOG:001.1a.2c'))]
        for hog, rng in cases:
            self.assertEqual(rng, self.db._hog_lex_range(hog))

    def test_fam_member(self):
        memb = self.db.member_of_fam(1)
        self.assertEqual(2, len(memb))

    def test_hog_members(self):
        cases = [('Eukaryota', 5), ('Fungi', 3), ('Taphrinomycotina',0)]
        for level, exp_member_cnt in cases:
            if exp_member_cnt == 0:
                with self.assertRaises(ValueError):
                    self.db.hog_members(12, level)
            else:
                self.assertEqual(len(self.db.hog_members(12, level)), exp_member_cnt)

    def test_hogids_at_level(self):
        cases = [[(2, 'Ascomycota'), numpy.array([b'HOG:0000002'])],
                 [(2, 'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)'),
                  numpy.array([b'HOG:0000002.2a', b'HOG:0000002.2b'])]]

        for case in cases:
            args, expected = case
            levels = self.db.get_subhogids_at_level(*args)
            self.assertTrue(numpy.array_equal(expected, levels),
                            'test of tes_hogids_at_level failed for {}: {}'.format(args, levels))


class TaxonomyTest(unittest.TestCase):
    tax_input = None
    maxDiff = None  # show complete strings if test fails

    @classmethod
    def setUpClass(cls):
        h5 = tables.open_file('/pub/projects/cbrg-oma-browser/Test.Jul2014/data/OmaServer.h5')
        cls.tax_input = h5.root.Taxonomy.read()
        h5.close()

    def setUp(self):
        self.tax = Taxonomy(self.tax_input)

    def test_parents(self):
        lin = [x['Name'].decode() for x in self.tax.get_parent_taxa(284811)]
        self.assertEqual(lin, ['Ashbya gossypii (strain ATCC 10895 / CBS 109.51 / FGSC 9923 / NRRL Y-1056)',
                               'Eremothecium', 'Saccharomycetaceae', 'Saccharomycetales',
                               'Saccharomycetes', 'Saccharomycotina', 'saccharomyceta',
                               'Ascomycota', 'Dikarya', 'Fungi', 'Opisthokonta', 'Eukaryota'])

    def test_newick(self):
        member = frozenset([self.tax._taxon_from_numeric(x)['Name']
                            for x in self.tax.tax_table['NCBITaxonId']])
        phylo = self.tax.get_induced_taxonomy(member, collapse=True)
        expected = '(((Ashbya gossypii (strain ATCC 10895 / CBS 109.51 / FGSC 9923 / NRRL Y-1056),Saccharomyces cerevisiae (strain ATCC 204508 / S288c))Saccharomycetaceae,Schizosaccharomyces pombe (strain 972 / ATCC 24843))Ascomycota,Plasmodium falciparum (isolate 3D7))Eukaryota'
        expected = expected.replace(' ', '_')
        self.assertEqual(expected, phylo.newick())

    def test_phylogeny(self):
        member = frozenset([self.tax._taxon_from_numeric(x)['Name']
                            for x in self.tax.tax_table['NCBITaxonId']])
        phylo = self.tax.get_induced_taxonomy(member, collapse=True)
        expected = {"id":2759, "name":"Eukaryota","children":[
                      {"id":4890, "name":"Ascomycota", "children":[
                        {"id":4893, "name": "Saccharomycetaceae", "children":[
                          {"id":284811, "name": "Ashbya gossypii (strain ATCC 10895 / CBS 109.51 / FGSC 9923 / NRRL Y-1056)"},
                          {"id":559292, "name": "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"}]},
                        {"id":284812, "name": "Schizosaccharomyces pombe (strain 972 / ATCC 24843)"}]},
                      {"id":36329, "name": "Plasmodium falciparum (isolate 3D7)"}]}
        self.assertEqual(expected, phylo.as_dict())
