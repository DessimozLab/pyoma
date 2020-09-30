# This Python file uses the following encoding: latin-1
from __future__ import unicode_literals

import tempfile
import unittest

try:
    import unittest.mock as mock
except ImportError:
    import mock as mock
import numpy
import tables
import io
import os
from pyoma.browser import convert as pyoma


@mock.patch.object(
    pyoma.UniProtAdditionalXRefImporter,
    "iter_xreftuples_for_up",
    autospec=True,
    return_value=[],
)
class XRefParsingTest(unittest.TestCase):
    def _create_genome_info(self):
        gs = numpy.zeros(4, dtype=tables.dtype_from_descr(pyoma.tablefmt.GenomeTable))
        gs["EntryOff"] = [0, 2, 3, 4]
        gs["TotEntries"] = [2, 1, 1, 1]
        gs["Release"] = [
            b"Ensembl",
            b"madeup",
            b"Ensembl Metazoa 22; CB4; 13-MAR-2014",
            b"Ensembl Metazoa 27; GCA_000002325.2; 2-JUN-2015",
        ]
        return gs

    def setUp(self):
        self.data = io.StringIO(
            """<E><ID>ENSG00000204640</ID><AC>ENSP00000366061; ENST00000376865</AC><DE>hypotetical protein</DE><GI>125233342</GI><UniProt/TrEMBL>P21522</UniProt/TrEMBL></E>
            <E><ID>BLA22; BLABLA22.Rep22</ID><AC>BLA22.1</AC><EC>3.2.2.-; 4.2.99.18</EC><EntrezGene>32244</EntrezGene><SMR>P21122; Q24S32</SMR><GO>GO:0006270@[[IDA,{'PMID:2167836'}],[IEA,{'GO_REF:002','GO_REF:020','OMA_Fun:001'}]]; GO:0006275@[[IEA,{'GO_REF:002','GO_REF:020','OMA_Fun:001'}]]</GO></E>
            <E><UniProt/TrEMBL>L8ECQ9_BACSU</UniProt/TrEMBL><SwissProt_AC>Q6CI62</SwissProt_AC><SwissProt>ASF1_YARLI</SwissProt><ID>FBgn0218776</ID><AC>FBpp0245919; FBtr0247427</AC><DE>β-hemoglobin</DE><GO></GO></E>
            <E><OS>CAEBR</OS><NR>1</NR><OG>0</OG><AC>CBG23988</AC><CHR>chrI</CHR><ID>CBG23988</ID><LOC>join(80..1057,2068..2664)</LOC><UniProt/TrEMBL>A8WJQ9_CAEBR</UniProt/TrEMBL><GO>GO:0016020@[[IEA,{'GO_REF:038'}]]; GO:0016021@[[IEA,{'GO_REF:038'}]]</GO><SEQ>A</SEQ></E>
            <E><OS>NASVI</OS><NR>11661</NR><OG>0</OG><AC>NV21158-PA; NV21158-RA</AC><CHR>GL340889</CHR><ID>NV21158</ID><LOC>join(113393..113590,114511..114696,114865..114916,115242..115474)</LOC><UniProt/TrEMBL>K7JG62_NASVI</UniProt/TrEMBL><ProtName>Uncharacterized protein</ProtName><SEQ>A</SEQ></E>"""
        )
        self.db_parser = pyoma.DarwinDbEntryParser()
        self.desc_manager = mock.Mock()
        self.go_manager = mock.Mock()
        self.approx_adder = pyoma.ApproximateXRefImporter()
        self.up_extra_adder = pyoma.UniProtAdditionalXRefImporter()
        self.xref_tab = mock.Mock()
        self.ec_tab = mock.Mock()
        self.gs = self._create_genome_info()
        self.importer = pyoma.XRefImporter(
            self.db_parser,
            self.gs,
            self.xref_tab,
            self.ec_tab,
            self.go_manager,
            self.desc_manager,
            self.approx_adder,
            self.up_extra_adder,
        )

    def test_standard_handler(self, mock_up):
        self.db_parser.parse_entrytags(self.data)
        enum = pyoma.tablefmt.XRefTable.columns.get("XRefSource").enum
        verif = pyoma.tablefmt.XRefTable.columns.get("Verification").enum
        self.assertIn((1, enum.GI, b"125233342", verif.exact), self.importer.xrefs)
        self.assertIn(
            (3, enum["UniProtKB/SwissProt"], b"ASF1_YARLI", verif.exact),
            self.importer.xrefs,
        )

    def test_multi_handler(self, mock_up):
        self.db_parser.parse_entrytags(self.data)
        enum = pyoma.tablefmt.XRefTable.columns.get("XRefSource").enum
        verif = pyoma.tablefmt.XRefTable.columns.get("Verification").enum
        self.assertIn(
            (2, enum["Swiss Model"], b"P21122", verif.unchecked), self.importer.xrefs
        )
        self.assertIn(
            (2, enum["Swiss Model"], b"Q24S32", verif.unchecked), self.importer.xrefs
        )

    def test_regex_of_ensembl_ids(self, mock_up):
        for case in ("ENSG00000162687", "ENSMUSP00000162687"):
            match = self.importer.ENS_RE.match(case)
            self.assertIsNotNone(match)

    def test_potential_flush_gets_called(self, mock_up):
        callback = mock.MagicMock(name="potential_flush")
        self.db_parser.add_post_entry_handler(callback)
        self.db_parser.parse_entrytags(self.data)
        self.assertEqual(callback.call_count, 5)

    def test_ensembgenomes_ids(self, mock_up):
        self.db_parser.parse_entrytags(self.data)
        enum = pyoma.tablefmt.XRefTable.columns.get("XRefSource").enum
        verif = pyoma.tablefmt.XRefTable.columns.get("Verification").enum
        self.assertIn(
            (4, enum.EnsemblGenomes, b"CBG23988", verif.exact), self.importer.xrefs
        )

    def test_uniprot_ids(self, mock_up):
        self.db_parser.parse_entrytags(self.data)
        enum = pyoma.tablefmt.XRefTable.columns.get("XRefSource").enum
        verif = pyoma.tablefmt.XRefTable.columns.get("Verification").enum
        self.assertIn(
            (3, enum["UniProtKB/TrEMBL"], b"L8ECQ9", verif.exact), self.importer.xrefs
        )
        self.assertIn(
            (3, enum["UniProtKB/TrEMBL"], b"Q6CI62", verif.exact), self.importer.xrefs
        )
        mock_up.assert_called_with(self.up_extra_adder, "K7JG62", 5)

    def test_go(self, mock_up):
        self.db_parser.parse_entrytags(self.data)
        self.assertEqual(len(self.go_manager.add_annotations.call_args_list), 2)

    def test_ec(self, mock_up):
        self.db_parser.parse_entrytags(self.data)
        self.assertIn((2, "3.2.2.-"), self.importer.ec)
        self.assertIn((2, "4.2.99.18"), self.importer.ec)

    def test_nasvi_missing_ac(self, mock_up):
        self.db_parser.parse_entrytags(self.data)
        self.importer.flush_buffers()
        enum = pyoma.tablefmt.XRefTable.columns.get("XRefSource").enum
        verif = pyoma.tablefmt.XRefTable.columns.get("Verification").enum
        args, kwargs = self.importer.xref_tab.append.call_args_list[0]
        buffer = args[0]
        self.assertIn(
            (5, enum.SourceAC, b"NV21158-PA", verif.exact), buffer,
        )
        self.assertIn(
            (5, enum.SourceAC, b"NV21158-RA", verif.exact), buffer,
        )
        self.assertIn(
            (5, enum.EnsemblGenomes, b"NV21158-RA", verif.exact), buffer,
        )

    def test_disambiguate(self, mock_up):
        self.db_parser.parse_entrytags(self.data)
        enum = pyoma.tablefmt.XRefTable.columns.get("XRefSource").enum
        verif = pyoma.tablefmt.XRefTable.columns.get("Verification").enum
        self.assertIn(
            (3, enum["FlyBase"], b"FBgn0218776", verif.unchecked), self.importer.xrefs
        )
        self.assertIn(
            (3, enum["FlyBase"], b"FBtr0247427", verif.unchecked), self.importer.xrefs
        )
        self.assertIn(
            (3, enum["SourceAC"], b"FBtr0247427", verif.exact), self.importer.xrefs
        )
        self.assertIn(
            (2, enum["SourceID"], b"BLABLA22", verif.exact), self.importer.xrefs
        )
        self.assertIn((2, enum["SourceID"], b"BLA22", verif.exact), self.importer.xrefs)
        self.assertIn(
            (1, enum["SourceID"], b"ENSG00000204640", verif.exact), self.importer.xrefs
        )
        self.assertIn(
            (1, enum["Ensembl Gene"], b"ENSG00000204640", verif.exact),
            self.importer.xrefs,
        )
        self.assertIn(
            (1, enum["Ensembl Protein"], b"ENSP00000366061", verif.exact),
            self.importer.xrefs,
        )
        self.assertIn(
            (1, enum["Ensembl Transcript"], b"ENST00000376865", verif.exact),
            self.importer.xrefs,
        )

    def test_descriptions_passed_to_description_manager(self, mock_up):
        self.db_parser.parse_entrytags(self.data)
        self.assertEqual(len(self.desc_manager.add_description.call_args_list), 2)
        self.assertEqual(
            (3, "β-hemoglobin"), self.desc_manager.add_description.call_args[0]
        )

    def test_remove_duplicated_xrefs(self, mock_up):
        ref = (1, 10, "test_id", "unchecked")
        res = (1, 10, b"test_id", 2)
        self.importer.FLUSH_SIZE = 20
        for i in range(self.importer.FLUSH_SIZE + 1):
            self.importer._add_to_xrefs(*ref)
        self.xref_tab.append.assert_not_called()
        self.importer.potential_flush()
        self.xref_tab.append.assert_called_once_with([res])

    def test_remove_near_duplicated_xrefs(self, mock_up):
        refs = [(1, 10, "test_id", "unchecked"), (1, 10, "test_id", "exact")]
        res = (1, 10, b"test_id", 2)
        self.importer.FLUSH_SIZE = 20
        for i in range(self.importer.FLUSH_SIZE + 1):
            for ref in refs:
                self.importer._add_to_xrefs(*ref)
        self.xref_tab.append.assert_not_called()
        self.importer.potential_flush()
        self.xref_tab.append.assert_called_once_with([res])


class DescriptionManagerTest(unittest.TestCase):
    def setUp(self):
        h5file = tables.open_file(
            "test.h5", "w", driver="H5FD_CORE", driver_core_backing_store=0
        )
        nr_rows = 3
        data = numpy.zeros(
            nr_rows, dtype=tables.dtype_from_descr(pyoma.tablefmt.ProteinTable)
        )
        data["EntryNr"] = numpy.arange(1, nr_rows + 1)
        h5file.create_table("/", "Entries", pyoma.tablefmt.ProteinTable, obj=data)
        self.h5 = h5file

    def tearDown(self):
        self.h5.close()

    def check_descriptions_one_per_entry(self, descs):
        with pyoma.DescriptionManager(self.h5, "/Entries", "/Descr") as dm:
            for i, desc in enumerate(descs):
                dm.add_description(i + 1, desc)
        exp_lens = numpy.array(list(len(z.encode("utf-8")) for z in descs), dtype="i4")
        exp_offs = numpy.cumsum(exp_lens) - exp_lens[0]
        entry_tab = self.h5.get_node("/Entries")
        numpy.array_equal(exp_lens, entry_tab.col("DescriptionLength"))
        numpy.array_equal(exp_offs, entry_tab.col("DescriptionOffset"))
        txt = self.h5.get_node("/Descr").read().tostring()
        self.assertEqual("".join(descs).encode("utf-8"), txt)

    def test_add_one_desc_per_entry(self):
        descs = ["hello", "some more lengthy test", "a"]
        self.check_descriptions_one_per_entry(descs)

    def test_add_with_unicode(self):
        descs = ["β-hemoglobin", "", "foobar"]
        self.check_descriptions_one_per_entry(descs)

    def test_severl_descr_per_entry(self):
        descs = [["blabla", "βgamma", "bar"], [""], ["xxx", "yyy"]]
        with pyoma.DescriptionManager(self.h5, "/Entries", "/Descr") as dm:
            for i, per_entry_set in enumerate(descs):
                for desc in per_entry_set:
                    dm.add_description(i + 1, desc)
        exp_txt = "".join(["; ".join(z) for z in descs])
        txt_in_db = self.h5.get_node("/Descr").read().tostring()
        self.assertEqual(exp_txt.encode("utf-8"), txt_in_db)


class TestGeneOntologyManager(pyoma.GeneOntologyManager):
    ontology_url = "file://" + os.path.dirname(__file__) + "/go-basic-tiny.obo"


class GeneOntologyManagerTest(unittest.TestCase):
    def setUp(self):
        self.h5file = tables.open_file(
            "test.h5", "w", driver="H5FD_CORE", driver_core_backing_store=0
        )

    def tearDown(self):
        self.h5file.close()

    def store_and_retrieve_data(self, enr, gos):
        with TestGeneOntologyManager(self.h5file, "/Annotation/GO", "/obo") as man:
            man.add_annotations(enr, gos)
        return self.h5file.get_node("/Annotation/GO").read()

    def test_add_none_as_data_raises_ValueError(self):
        with self.assertRaises(ValueError):
            self.store_and_retrieve_data(None, None)

    def test_basic_anno(self):
        res = self.store_and_retrieve_data(
            1,
            "GO:0007610@[[IDA,{'PMID:2167836'}],[IEA,{'GO_REF:002','GO_REF:020','OMA_Fun:001'}]]; GO:0019954@[[IEA,{'GO_REF:002','GO_REF:020','OMA_Fun:001'}]]",
        )
        self.assertEqual(7, len(res))
        numpy.testing.assert_equal(1, res["EntryNr"])
        self.assertEqual(7610, res[0]["TermNr"])
        self.assertEqual(b"IDA", res[0]["Evidence"])
        self.assertEqual(b"PMID:2167836", res[0]["Reference"])

    def test_obo_version_set(self):
        self.store_and_retrieve_data(2, "")
        go_node = self.h5file.get_node("/obo")
        vers = go_node.attrs["ontology_release"]
        self.assertEqual("dummy-test/2017-04-18", vers)

    def test_add_obsolete_term_is_skipped(self):
        import sys

        if sys.version_info >= (3, 4):
            with self.assertLogs("pyoma", level="WARNING") as log:
                res = self.store_and_retrieve_data(
                    1, "GO:0003822@[[IEA,{OMA_Fun:001}]]"
                )
                self.assertIn(
                    "invalid GO term for entry 1: GO:0003822", ";".join(log.output)
                )
        else:
            res = self.store_and_retrieve_data(1, "GO:0003822@[[IEA,{OMA_Fun:001}]]")
        self.assertEqual(0, len(res))

    def test_go_obo_version(self):
        self.store_and_retrieve_data(1, "")
        obo = self.h5file.get_node("/obo")
        self.assertEqual("dummy-test/2017-04-18", obo._f_getattr("ontology_release"))


class ApproxXRefImporterTest(unittest.TestCase):
    def setUp(self):
        with tempfile.NamedTemporaryFile(delete=False) as fh:
            self.approx_fname = fh.name
            fh.write(
                """1315\tGeneName\tFca\tmodified
1315\tProtName\tFCA protein\tmodified
1315\tUniProtKB/TrEMBL\tQ6XJS3\tmodified
1406\tORFNames\tAALP_AA5G153500\tmodified
1229\tORFNames\tL914_11153\tmodified""".encode(
                    "utf-8"
                )
            )
        self.approx = pyoma.ApproximateXRefImporter(self.approx_fname)

    def test_non_existing(self):
        self.assertEqual(len(list(self.approx.iter_approx_xrefs_for(10))), 0)

    def test_exitsting(self):
        res = list(self.approx.iter_approx_xrefs_for(1229))
        self.assertEqual(1, len(res), res)
        self.assertEqual(res[0][0], 1229)

    def test_non_sorted_access_raises(self):
        list(self.approx.iter_approx_xrefs_for(100))
        with self.assertRaises(pyoma.InvarianceException):
            list(self.approx.iter_approx_xrefs_for(90))

    def test_return_all_for_large_range(self):
        res = []
        for enr in range(2000):
            res.extend(self.approx.iter_approx_xrefs_for(enr))
        self.assertEqual(len(self.approx.xrefs), len(res))

    def tearDown(self):
        os.remove(self.approx_fname)
