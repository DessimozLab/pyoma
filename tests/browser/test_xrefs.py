import unittest
import unittest.mock

import numpy
import tables
import io
from pyoma.browser import convert as pyoma


class XRefParsingTest(unittest.TestCase):
    def setUp(self):
        data = io.StringIO(
            """<E><ID>ENSG00000204640</ID><AC>ENSP00000366061; ENST00000376865</AC><DE>hypotetical protein</DE><GI>125233342</GI><UniProt/TrEMBL>P21522</UniProt/TrEMBL></E>
            <E><ID>BLA22; BLABLA22.Rep22</ID><AC>BLA22.1</AC><EntrezGene>32244</EntrezGene><PMP>P21122; Q24S32</PMP><GO>GO:0006270@[[IDA,{'PMID:2167836'}],[IEA,{'GO_REF:002','GO_REF:020','OMA_Fun:001'}]]; GO:0006275@[[IEA,{'GO_REF:002','GO_REF:020','OMA_Fun:001'}]]</GO></E>
            <E><UniProt/TrEMBL>L8ECQ9_BACSU</UniProt/TrEMBL><SwissProt_AC>Q6CI62</SwissProt_AC><SwissProt>ASF1_YARLI</SwissProt><ID>FBgn0218776</ID><AC>FBpp0245919; FBtr0247427</AC><DE>β-hemoglobin</DE><GO></GO></E>""")
        self.db_parser = pyoma.IndexedServerParser(data)
        self.desc_manager = unittest.mock.Mock()
        self.importer = pyoma.XRefImporter(self.db_parser, None, None, self.desc_manager)

    def test_standard_handler(self):
        self.db_parser.parse_entrytags()
        enum = pyoma.tablefmt.XRefTable.columns.get('XRefSource').enum
        self.assertIn((1, enum.GI, b'125233342',), self.importer.xrefs)
        self.assertIn((3, enum['UniProtKB/SwissProt'], b'ASF1_YARLI',), self.importer.xrefs)

    def test_multi_handler(self):
        self.db_parser.parse_entrytags()
        enum = pyoma.tablefmt.XRefTable.columns.get('XRefSource').enum
        self.assertIn((2, enum.PMP, b'P21122',), self.importer.xrefs)
        self.assertIn((2, enum.PMP, b'Q24S32',), self.importer.xrefs)

    def test_regex_of_ensembl_ids(self):
        for case in ('ENSG00000162687', 'ENSMUSP00000162687'):
            match = self.importer.ENS_RE.match(case)
            self.assertIsNotNone(match)

    def test_uniprot_ids(self):
        self.db_parser.parse_entrytags()
        enum = pyoma.tablefmt.XRefTable.columns.get('XRefSource').enum
        self.assertIn((3, enum['UniProtKB/TrEMBL'], b'L8ECQ9',), self.importer.xrefs)
        self.assertIn((3, enum['UniProtKB/TrEMBL'], b'Q6CI62',), self.importer.xrefs)

    def test_go(self):
        self.db_parser.parse_entrytags()
        self.assertIn((2, 6270, 'IDA', b'PMID:2167836'), self.importer.go)
        self.assertIn((2, 6275, 'IEA', b'OMA_Fun:001'), self.importer.go)

    def test_disambiguate(self):
        self.db_parser.parse_entrytags()
        enum = pyoma.tablefmt.XRefTable.columns.get('XRefSource').enum
        self.assertIn((3, enum['FlyBase'], b'FBgn0218776'), self.importer.xrefs)
        self.assertIn((3, enum['FlyBase'], b'FBtr0247427'), self.importer.xrefs)
        self.assertIn((3, enum['SourceAC'], b'FBtr0247427'), self.importer.xrefs)
        self.assertIn((2, enum['SourceID'], b'BLABLA22'), self.importer.xrefs)
        self.assertIn((2, enum['SourceID'], b'BLA22'), self.importer.xrefs)
        self.assertIn((1, enum['SourceID'], b'ENSG00000204640'), self.importer.xrefs)
        self.assertIn((1, enum['Ensembl Gene'], b'ENSG00000204640'), self.importer.xrefs)
        self.assertIn((1, enum['Ensembl Protein'], b'ENSP00000366061'), self.importer.xrefs)
        self.assertIn((1, enum['Ensembl Transcript'], b'ENST00000376865'), self.importer.xrefs)

    def test_descriptions_passed_to_description_manager(self):
        self.db_parser.parse_entrytags()
        self.assertEqual(len(self.desc_manager.add_description.call_args_list), 2)
        self.assertEqual((3,u'β-hemoglobin'), self.desc_manager.add_description.call_args[0])


class DescriptionManagerTest(unittest.TestCase):
    def setUp(self):
        h5file = tables.open_file("test.h5", "w", driver="H5FD_CORE",
                                  driver_core_backing_store=0)
        nr_rows = 3
        data = numpy.zeros(nr_rows, dtype=tables.dtype_from_descr(pyoma.tablefmt.ProteinTable))
        data['EntryNr'] = numpy.arange(1,nr_rows + 1)
        h5file.create_table('/', 'Entries', pyoma.tablefmt.ProteinTable, obj=data)
        self.h5 = h5file

    def tearDown(self):
        self.h5.close()

    def check_descriptions_one_per_entry(self, descs):
        with pyoma.DescriptionManager(self.h5, '/Entries', '/Descr') as dm:
            for i, desc in enumerate(descs):
                dm.add_description(i + 1, desc)
        exp_lens = numpy.array(list(len(z.encode('utf-8')) for z in descs), dtype='i4')
        exp_offs = numpy.cumsum(exp_lens) - exp_lens[0]
        entry_tab = self.h5.get_node('/Entries')
        numpy.array_equal(exp_lens, entry_tab.col('DescriptionLength'))
        numpy.array_equal(exp_offs, entry_tab.col('DescriptionOffset'))
        txt = self.h5.get_node('/Descr').read().tostring()
        self.assertEqual("".join(descs).encode('utf-8'), txt)

    def test_add_one_desc_per_entry(self):
        descs = ["hello", "some more lengthy test", "a"]
        self.check_descriptions_one_per_entry(descs)

    def test_add_with_unicode(self):
        descs = ["β-hemoglobin", "", "foobar"]
        self.check_descriptions_one_per_entry(descs)

    def test_severl_descr_per_entry(self):
        descs = [["blabla", "βgamma", "bar"], [""], ["xxx", "yyy"]]
        with pyoma.DescriptionManager(self.h5, '/Entries', '/Descr') as dm:
            for i, per_entry_set in enumerate(descs):
                for desc in per_entry_set:
                    dm.add_description(i + 1, desc)
        exp_txt = "".join(["; ".join(z) for z in descs])
        txt_in_db = self.h5.get_node('/Descr').read().tostring()
        self.assertEqual(exp_txt.encode('utf-8'), txt_in_db)