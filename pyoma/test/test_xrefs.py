import unittest
try:
    import io
except ImportError:
    import stringIO as io
from .. import convert as pyoma


class XRefParsingTest(unittest.TestCase):

    def setUp(self):
        data = io.StringIO(
            """<E><ID>ENSG00000204640</ID><AC>ENSP00000366061; ENST00000376865</AC><DE>hypotetical protein</DE><GI>125233342</GI><UniProt/TrEMBL>P21522</UniProt/TrEMBL></E>
            <E><ID>BLA22; BLABLA22.Rep22</ID><AC>BLA22.1</AC><EntrezGene>32244</EntrezGene><PMP>P21122; Q24S32</PMP><GO>GO:0006270@[[IDA,{'PMID:2167836'}],[IEA,{'GO_REF:002','GO_REF:020','OMA_Fun:001'}]]; GO:0006275@[[IEA,{'GO_REF:002','GO_REF:020','OMA_Fun:001'}]]</GO></E>
            <E><UniProt/TrEMBL>L8ECQ9_BACSU</UniProt/TrEMBL><SwissProt_AC>Q6CI62</SwissProt_AC><SwissProt>ASF1_YARLI</SwissProt><ID>FBgn0218776</ID><AC>FBpp0245919; FBtr0247427</AC><DE>β-hemoglobin</DE><GO></GO></E>""")
        self.db_parser = pyoma.IndexedServerParser(data)
        self.importer = pyoma.XRefImporter(self.db_parser)


    def test_standard_handler(self):
        self.db_parser.parse_entrytags()
        enum = pyoma.XRefTable.columns.get('XRefSource').enum
        self.assertIn((1, enum.GI, b'125233342',), self.importer.xrefs)
        self.assertIn((1, enum.Description, b'hypotetical protein',), self.importer.xrefs)
        self.assertIn((3, enum['UniProtKB/SwissProt'], b'ASF1_YARLI',), self.importer.xrefs)

    
    def test_multi_handler(self):
        self.db_parser.parse_entrytags()
        enum = pyoma.XRefTable.columns.get('XRefSource').enum
        self.assertIn((2, enum.PMP, b'P21122',), self.importer.xrefs)
        self.assertIn((2, enum.PMP, b'Q24S32',), self.importer.xrefs)

    def test_regex_of_ensembl_ids(self):
        for case in ('ENSG00000162687','ENSMUSP00000162687'):
            match = self.importer.ENS_RE.match(case)
            self.assertIsNotNone(match)

    def test_uniprot_ids(self):
        self.db_parser.parse_entrytags()
        enum = pyoma.XRefTable.columns.get('XRefSource').enum
        self.assertIn((3, enum['UniProtKB/TrEMBL'], b'L8ECQ9',), self.importer.xrefs)
        self.assertIn((3, enum['UniProtKB/TrEMBL'], b'Q6CI62',), self.importer.xrefs)
    
    def test_go(self):
        self.db_parser.parse_entrytags()
        self.assertIn((2, 6270, 'IDA', b'PMID:2167836'), self.importer.go)
        self.assertIn((2, 6275, 'IEA', b'OMA_Fun:001'), self.importer.go)

    def test_disambiguate(self):
        self.db_parser.parse_entrytags()
        enum = pyoma.XRefTable.columns.get('XRefSource').enum
        self.assertIn((3,enum['FlyBase'], b'FBgn0218776'), self.importer.xrefs)
        self.assertIn((3,enum['FlyBase'], b'FBtr0247427'), self.importer.xrefs)
        self.assertIn((3,enum['SourceAC'], b'FBtr0247427'), self.importer.xrefs)
        self.assertIn((2,enum['SourceID'], b'BLABLA22'), self.importer.xrefs) 
        self.assertIn((2,enum['SourceID'], b'BLA22'), self.importer.xrefs)
        self.assertIn((1,enum['SourceID'], b'ENSG00000204640'), self.importer.xrefs)
        self.assertIn((1,enum['Ensembl Gene'], b'ENSG00000204640'), self.importer.xrefs)
        self.assertIn((1,enum['Ensembl Protein'], b'ENSP00000366061'), self.importer.xrefs)
        self.assertIn((1,enum['Ensembl Transcript'], b'ENST00000376865'), self.importer.xrefs)

    def test_unicode_encoding(self):
        self.db_parser.parse_entrytags()
        enum = pyoma.XRefTable.columns.get('XRefSource').enum
        self.assertIn( (3,enum.Description, 'β-hemoglobin'.encode('utf-8')), self.importer.xrefs)

