from __future__ import division, print_function, unicode_literals

import random
import types
import unittest

import numpy.testing

try:
    import unittest.mock as mock
except ImportError:
    import mock
from pyoma.browser.db import *
from pyoma.browser import tablefmt
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


class TestHelperFunctions(unittest.TestCase):
    def test_counter(self):
        self.assertEqual(0, count_elements([]))
        self.assertEqual(3, count_elements("abc"))
        recarray = numpy.zeros(2, dtype=[("A", "i4"), ("B", "f8")])
        self.assertEqual(2, count_elements(recarray))


def find_path_to_test_db(dbfn="TestDb.h5"):
    """We try to load the dbfn first from the same directory, and afterwards from
    the path given by the PYOMA_DB_PATH environment variable.

    :returns: path to database
    :rtype: str
    :raises IOError: if db does not exist."""

    def is_hdf5_file(fname):
        with open(fname, "rb") as fh:
            return fh.read(4) == b"\x89HDF"

    path1 = os.path.join(os.path.dirname(__file__), dbfn)
    if os.path.isfile(path1) and is_hdf5_file(path1):
        return path1
    path2 = os.path.abspath(os.path.join(os.getenv("PYOMA_DB_PATH", "./"), dbfn))
    if os.path.isfile(path2) and is_hdf5_file(path2):
        return path2
    else:
        raise IOError("cannot access {}. (Tried {} and {})".format(dbfn, path1, path2))


class TestWithDbInstance(unittest.TestCase):
    db = None
    db_file_name = "TestDb.h5"

    @classmethod
    def setUpClass(cls):
        path = find_path_to_test_db(cls.db_file_name)
        logger.info("Loading {} for DatabaseTests".format(path))
        cls.db = Database(path)

    @classmethod
    def tearDownClass(cls):
        cls.db.close()


class DatabaseTests(TestWithDbInstance):
    def test_get_vpairs_of_entry_with_orthologs(self):
        for entry_nr, exp_vps_cnt in [(12, 3), (1, 0), (4, 1)]:
            vps = self.db.get_vpairs(entry_nr)
            self.assertTrue(isinstance(vps, numpy.ndarray))
            self.assertEqual(exp_vps_cnt, len(vps))
            self.assertEqual(
                sorted(["EntryNr1", "EntryNr2", "RelType", "Distance", "Score"]),
                sorted(vps.dtype.fields.keys()),
            )

    def test_get_vpairs_between_species_is_same_as_loading_all_vpairs(self):
        g1, g2 = "YEAST", "SCHPO"
        t0 = time.time()
        rng1 = self.db.id_mapper["OMA"].genome_range(g1)
        rng2 = self.db.id_mapper["OMA"].genome_range(g2)
        expected = []
        for en in range(rng1[0], rng1[1] + 1):
            for vp in self.db.get_vpairs(en):
                if rng2[0] <= vp["EntryNr2"] <= rng2[1]:
                    expected.append(vp)
        expected = numpy.fromiter(expected, dtype=expected[0].dtype)
        t1 = time.time()
        actual = self.db.get_vpairs_between_species_pair(g1, g2)
        t2 = time.time()
        numpy.testing.assert_equal(expected, actual)
        self.assertLess(10 * (t2 - t1), t1 - t0, "expect at least 10x speedup")

    def test_get_hoginducedpairs_between_species_is_same_as_loading_all_hoginduced(
        self,
    ):
        g1, g2 = "YEAST", "SCHPO"
        t0 = time.time()
        rng1 = self.db.id_mapper["OMA"].genome_range(g1)
        rng2 = self.db.id_mapper["OMA"].genome_range(g2)
        expected = []
        for en in range(rng1[0], rng1[1] + 1):
            for vp in self.db.get_hog_induced_pairwise_orthologs(en):
                if rng2[0] <= vp["EntryNr"] <= rng2[1]:
                    expected.append((en, int(vp["EntryNr"])))
        expected = numpy.fromiter(
            expected, dtype=[("EntryNr1", "u4"), ("EntryNr2", "u4")]
        )
        t1 = time.time()
        actual = self.db.get_hog_induced_orthologs_between_genome_pair(g1, g2)
        actual = actual[["EntryNr1", "EntryNr2"]].copy()
        actual.sort()
        t2 = time.time()
        numpy.testing.assert_equal(expected, actual)
        self.assertLess(5 * (t2 - t1), t1 - t0, "expect at least 5x speedup")

    def test_neighborhood_close_to_boundary(self):
        query, window = 3, 7
        neighbors, idx = self.db.neighbour_genes(query, window)
        self.assertEqual(query - 1, idx)
        self.assertEqual(neighbors["EntryNr"][idx], query)
        expected_entry_nrs = numpy.arange(1, query + window + 1, dtype="i4")
        self.assertTrue(numpy.array_equal(expected_entry_nrs, neighbors["EntryNr"]))

    def test_count_rows_of_index_column_with_value(self):
        entrynr = 12
        vps = self.db.get_vpairs(entrynr)
        vp_tab = self.db._get_vptab(entrynr)
        self.assertEqual(
            len(vps), count_rows_of_index_column_with_value(vp_tab, "EntryNr1", entrynr)
        )

        # check that it raises KeyError if operated on a column that is not indexed
        with self.assertRaises(
            KeyError, msg="Operating on non-indexed column should raise KeyError"
        ):
            count_rows_of_index_column_with_value(vp_tab, "EntryNr2", entrynr)

    def test_hog_family(self):
        entry = numpy.zeros(1, dtype=tables.dtype_from_descr(tablefmt.ProteinTable))
        entry["OmaHOG"] = b""
        with self.assertRaises(Singleton):
            self.db.hog_family(entry[0])
        entry["OmaHOG"] = b"HOG:0000523"
        self.assertEqual(523, self.db.hog_family(entry[0]))

    def test_orthoxml(self):
        xml = self.db.get_orthoxml(33).decode()
        # we simply check that orthoxml starts with <?xml and ends with an orthoxml tag
        self.assertTrue(xml.startswith("<?xml "))
        self.assertTrue(xml.endswith("</orthoXML>\n"))

    def test_hog_lex_range(self):
        cases = [
            (b"HOG:001", (b"HOG:001", b"HOG:002")),
            (b"HOG:001.1a.2b", (b"HOG:001.1a.2b", b"HOG:001.1a.2c")),
        ]
        for hog, rng in cases:
            self.assertEqual(rng, self.db._hog_lex_range(hog))

    def test_fam_member(self):
        memb = self.db.member_of_fam(1)
        self.assertEqual(2, len(memb))

    def test_hog_members(self):
        cases = [("Eukaryota", 5), ("Fungi", 3), ("Taphrinomycotina", 0)]
        for level, exp_member_cnt in cases:
            if exp_member_cnt == 0:
                with self.assertRaises(ValueError):
                    self.db.hog_members(12, level)
            else:
                self.assertEqual(len(self.db.hog_members(12, level)), exp_member_cnt)

    def test_hogids_at_level(self):
        cases = [
            [(2, "Ascomycota"), numpy.array([b"HOG:0000002"])],
            [
                (2, "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"),
                numpy.array([b"HOG:0000002.2a", b"HOG:0000002.2b"]),
            ],
        ]

        for case in cases:
            args, expected = case
            levels = self.db.get_subhogids_at_level(*args)
            self.assertTrue(
                numpy.array_equal(expected, levels),
                "test of tes_hogids_at_level failed for {}: {}".format(args, levels),
            )

    def test_get_hog_without_level(self):
        for case in ("HOG:0000002.2b", "HOG:0000001"):
            with self.subTest(case=case):
                hog = self.db.get_hog(case)
                self.assertTrue(hog["IsRoot"])
                self.assertEqual(case, hog["ID"].decode())

    def test_get_parent_hogs(self):
        for hog_id, level in (
            ("HOG:0000474.1d", "Ascomycota"),
            ("HOG:0000165", "Fungi"),
        ):
            with self.subTest(hog_id=hog_id, level=level):
                hogs = self.db.get_parent_hogs(hog_id, level)
                try:
                    k = hog_id.index(".")
                    base_id = hog_id[:k]
                except ValueError:
                    base_id = hog_id
                self.assertEqual(base_id, hogs[0].hog_id)
                self.assertTrue(hogs[0].is_root)
                self.assertEqual(hog_id, hogs[-1].hog_id)
                self.assertEqual(level, hogs[-1].level)

    def test_get_subhogs(self):
        for hog_id, level in (
            ("HOG:0000474.1d", "Ascomycota"),
            ("HOG:0000165", "Fungi"),
        ):
            with self.subTest(hog_id=hog_id, level=level):
                hogs = list(self.db.get_subhogs(hog_id, level))
                self.assertEqual(hogs[0].hog_id, hog_id)
                self.assertNotIn(level, [h.level for h in hogs])

    def test_get_subhogs_with_subids(self):
        for hog_id, level, exp_subhogid in (
            ("HOG:0000474", "Eukaryota", "HOG:0000474.1d.2c"),
            ("HOG:0000165", "Fungi", "HOG:0000165.1a"),
        ):
            with self.subTest(hog_id=hog_id, level=level):
                hogs = list(
                    self.db.get_subhogs(hog_id, level=level, include_subids=True)
                )
                self.assertIn(exp_subhogid, [h.hog_id for h in hogs])
                self.assertNotIn(level, [h.level for h in hogs])

    def test_get_hog_with_level(self):
        for hog_id, level in (
            ("HOG:0000002.2b", "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"),
            ("HOG:0000005.1b", "Saccharomycetaceae"),
            ("HOG:0000005", "Fungi"),
        ):
            with self.subTest(hog_id=hog_id, level=level):
                hog = self.db.get_hog(hog_id, level=level)
                self.assertEqual(hog_id, hog["ID"].decode())
                self.assertEqual(level, hog["Level"].decode())

    def test_member_of_hog_id(self):
        cases = [
            [("HOG:0000082.1b", None), 2],
            [("HOG:0000082.1b", "Saccharomycetaceae"), 2],
            [
                (
                    "HOG:0000082.1b",
                    "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)",
                ),
                1,
            ],
            [("HOG:0000082.1a", "Mammalia"), 0],
        ]
        for args, expected_len in cases:
            members = self.db.member_of_hog_id(*args)
            self.assertEqual(len(members), expected_len)

    def test_sorted_genomes(self):
        for root in ("YEAST", "ASHGO"):
            order = self.db.id_mapper["OMA"].species_ordering(root)
            self.assertEqual(
                order[root],
                0,
                "{} should be first genome, but comes at {}".format(root, order[root]),
            )

    def test_main_isoform_from_species_id(self):
        query = "PLAF7"
        rng = self.db.id_mapper["OMA"].genome_range(query)
        mains = self.db.main_isoforms(query)
        self.assertGreaterEqual(rng[1] - rng[0] + 1, len(mains))

    def test_valid_protein_seq_check(self):
        for c in (b"ADRIAN", "adrian", b"XaMGAt"):
            self.assertTrue(
                self.db.seq_search.contains_only_valid_chars(c),
                "'{}' should be a valid protein seq".format(c),
            )
        for c in (b"XBra", "T5HSE", "xxnnebqpynn", " "):
            self.assertFalse(
                self.db.seq_search.contains_only_valid_chars(c),
                "'{}' should be an invalid protein seq".format(c),
            )

    def test_exact_search(self):
        # Test for 10 random
        for _ in range(10):
            i = random.randint(0, len(self.db.db.root.Protein.Entries))
            enr = i + 1
            s = self.db.get_sequence(enr)
            self.assertTrue(
                (enr in set(self.db.seq_search.search(s, is_sanitised=True)[1])),
                "exact search for entry {} failed.".format(i),
            )

    def get_random_subsequence(self, minlen=10):
        i = random.randint(0, len(self.db.db.root.Protein.Entries))
        elen = self.db.db.root.Protein.Entries[i]["SeqBufferLength"] - 1
        enr = i + 1

        ii = random.randint(0, elen - minlen)
        jj = random.randint(ii + minlen, elen)

        s = self.db.get_sequence(enr)[ii:jj]
        return s, enr, ii, jj

    def test_approx_search(self):
        # Test for random subsequence of 10 random sequences.
        for _ in range(10):
            s, enr, start_idx, end_idx = self.get_random_subsequence()
            approx_search_results = self.db.seq_search.approx_search(
                s, is_sanitised=True
            )
            self.assertIn(
                enr,
                {z[0] for z in approx_search_results},
                "approx search for entry {}[{}:{}] failed.".format(
                    enr - 1, start_idx, end_idx
                ),
            )

    def test_specific_approx_search_that_failed_on_jenkins(self):
        # TODO fix this problem!
        enrs = [15885, 16452]
        ranges = [(39, 376), (55, 140)]
        for enr, (start_idx, end_idx) in zip(enrs, ranges):
            seq = self.db.get_sequence(enr)[start_idx:end_idx]
            approx_search_results = self.db.seq_search.approx_search(
                seq, is_sanitised=True
            )
            enrs_with_approx_match = {z[0] for z in approx_search_results}
            self.assertIn(enr, enrs_with_approx_match)

    def test_map_to_hog(self):
        hog_mapper = ClosestSeqMapper(self.db)
        seqs = []
        entry_nrs = []
        for _ in range(10):
            s, enr, start_idx, end_idx = self.get_random_subsequence()
            entry_nrs.append(enr)
            seqs.append(
                SeqRecord(
                    Seq(s.decode()),
                    id=str(enr),
                    annotations={"molecule_type": "protein"},
                )
            )
        res = list(
            filter(
                lambda x: x.query == str(x.closest_entry_nr),
                hog_mapper.imap_sequences(seqs),
            )
        )
        self.assertEqual(10, len(res))
        for case in res:
            self.assertTrue(0 < case.distance < 1)

    def test_oma_group_from_numeric_id(self):
        group_id = 5
        grp = self.db.oma_group_members(group_id)
        self.assertEqual(4, len(grp))
        for e in grp:
            self.assertEqual(group_id, e["OmaGroup"])

    def test_fingerprint(self):
        fingerprint, grp_nr = "ESRTELL", 2617
        grp = self.db.oma_group_members(fingerprint)
        self.assertLessEqual(2, len(grp))
        for e in grp:
            self.assertEqual(grp_nr, e["OmaGroup"])

    def test_exon_structure(self):
        query = 14677  # Q8I237
        exons = self.db.get_exons(query)
        self.assertEqual(3, len(exons))

    def test_go_term_search(self):
        query = "GO:0004575"
        nrs = self.db.entrynrs_with_go_annotation(query, evidence="IDA")
        self.assertGreaterEqual(
            len(nrs), 1, "query GO term is known to occure at least in MAL32_YEAST"
        )
        for enr in nrs:
            self.assertIn(4575, self.db.get_gene_ontology_annotations(enr)["TermNr"])

    def test_go_ic_values(self):
        terms = ["GO:003674", "GO:008150"]
        for term in terms:
            with self.subTest(case=term):
                gof = self.db.freq_aware_gene_ontology
                ic = gof.ic(term)
                self.assertAlmostEqual(0, ic)

    def test_go_lin_same_term_is_1(self):
        self.assertEqual(1, self.db.freq_aware_gene_ontology.lin_similarity(55, 55))

    def test_go_semantic_sim_of_same_term(self):
        self.assertGreater(
            1, self.db.freq_aware_gene_ontology.semantic_similarity(55, 55)
        )

    def test_mindepth_of_go(self):
        examples = [(8150, 0), (168, 7), (5689, 5)]
        for term, exp_depth in examples:
            with self.subTest(case=term):
                self.assertEqual(
                    exp_depth, self.db.gene_ontology.ensure_term(term).min_depth
                )

    def test_induced_pairwise_orthologs(self):
        query = "YEAST3523"
        query_entry = self.db.ensure_entry(self.db.id_resolver.resolve(query))
        orthologs = self.db.get_hog_induced_pairwise_orthologs(query_entry)
        self.assertEqual(3, len(orthologs))
        self.assertCountEqual([b"1:1", b"1:1", b"m:1"], orthologs["RelType"])

    def test_induced_pairwise_paralogs(self):
        query = "YEAST12"
        query_entry = self.db.ensure_entry(self.db.id_resolver.resolve(query))
        orthologs = self.db.get_hog_induced_pairwise_paralogs(query_entry)
        self.assertEqual(1, len(orthologs))
        self.assertEqual(
            b"Saccharomyces cerevisiae (strain ATCC 204508 / S288c)",
            orthologs["DivergenceLevel"],
        )

    def test_ancestral_synteny_graph(self):
        query_hog = "HOG:0000005"
        query_level = "Fungi"
        graph = self.db.get_syntentic_hogs(query_level, query_hog, steps=2)
        self.assertIn(query_hog, graph)

    def test_ancestral_syteny_of_genome(self):
        query_level = "Fungi"
        graph = self.db.get_syntenic_hogs(level=query_level)
        self.assertIn("HOG:0000005", graph)
        self.assertEqual(len(graph.nodes), self.db.count_hogs_at_level(query_level))

    def test_extant_syntey_of_genome(self):
        query_level = "YEAST"
        graph = self.db.get_extant_synteny_graph(query_level)
        self.assertIn("YEAST00012", graph)


class XRefDatabaseMock(Database):
    def __init__(self, name=None):
        if name is None:
            name = "xref.h5"
        f = tables.open_file(name, "w", driver="H5FD_CORE", driver_core_backing_store=0)
        entries = numpy.zeros(5, tables.dtype_from_descr(tablefmt.ProteinTable))
        entries["EntryNr"] = numpy.arange(1, 6)
        t = f.create_table("/Protein", "Entries", obj=entries, createparents=True)
        t.colinstances["EntryNr"].create_csindex()
        xref = numpy.zeros(10, tables.dtype_from_descr(tablefmt.XRefTable))
        xref["EntryNr"] = numpy.arange(1, 6, 0.5).astype(numpy.int32)
        xref["XRefSource"] = numpy.tile([0, 20], 5)
        xref["XRefId"] = ["XA{:03}g1.4".format(i) for i in range(10)]
        xref["Verification"] = tuple(itertools.islice(itertools.cycle([0, 2, 4]), 10))
        x = f.create_table("/", "XRef", tablefmt.XRefTable, obj=xref)
        x.colinstances["EntryNr"].create_csindex()
        x.colinstances["XRefId"].create_csindex()
        f.create_group("/", "XRef_Index")
        for n in ("suffix", "buffer", "offset"):
            f.create_carray("/XRef_Index", n, obj=numpy.ones((5,), "i4"))
        self.db = f
        self._on_close_notify = []
        self._close_fh = True


class XRefIdMapperTest(unittest.TestCase):
    @classmethod
    def setUp(self):
        patch_db = XRefDatabaseMock()
        self.xrefmapper = XrefIdMapper(patch_db)

    def tearDown(self):
        self.xrefmapper._db.close()

    def test_multiple_xrefs_per_entry(self):
        xref_e1 = self.xrefmapper.map_entry_nr(1)
        self.assertEqual(len(xref_e1), 2)

    def test_map_many_entries(self):
        all_mapped = self.xrefmapper.map_many_entry_nrs(numpy.arange(1, 4))
        expected_len = (
            4 if isinstance(self.xrefmapper, XRefNoApproximateIdMapper) else 6
        )
        self.assertEqual((expected_len,), all_mapped.shape)
        self.assertEqual(self.xrefmapper.xref_tab.dtype, all_mapped.dtype)

    def test_entry_nr_range(self):
        all_mapped = self.xrefmapper.map_many_entry_nrs(numpy.arange(1, 4))
        range_mapped = self.xrefmapper.map_entry_nr_range(1, 4)
        numpy.testing.assert_equal(all_mapped, range_mapped)

    def test_map_entry_iterator(self):
        it = self.xrefmapper.iter_xrefs_for_entry_nr(1)
        self.assertTrue(isinstance(it, types.GeneratorType), "not a generator")
        exp_xrefs = ["XA000g1.4", "XA001g1.4"]
        for dic in it:
            self.assertIn(dic["xref"], exp_xrefs)

    def test_search_of_modified_xref(self):
        xref = "XA002g1.4"
        res = self.xrefmapper.search_xref(xref)
        self.assertEqual(2, res["EntryNr"])

    def test_verification_is_returned(self):
        res = self.xrefmapper.map_entry_nr(2)
        if isinstance(self.xrefmapper, XRefNoApproximateIdMapper):
            expected = ["exact"]
        else:
            expected = ["modified", "exact"]
        self.assertEqual(expected, [z["seq_match"] for z in res])


class XRefNoApproxIdMapperTest(XRefIdMapperTest):
    def setUp(self):
        patch_db = XRefDatabaseMock()
        self.xrefmapper = XRefNoApproximateIdMapper(patch_db)

    def test_modified_xref_not_returned_in_map(self):
        res = self.xrefmapper.map_entry_nr(2)
        xrefs = [x["xref"] for x in res]
        self.assertNotIn("XA002g1.4", xrefs)
        self.assertIn("XA003g1.4", xrefs)


class IdResolverTests(unittest.TestCase):
    @classmethod
    def setUp(self):
        patch_db = XRefDatabaseMock("test_idresolver.h5")
        self.xrefmapper = XrefIdMapper(patch_db)
        self.id_resolver = IDResolver(patch_db)
        patch_db.id_mapper = {"XRef": self.xrefmapper}

    def tearDown(self):
        self.xrefmapper._db.close()

    def test_search_of_modified_xref(self):
        xref = "XA002g1.4"
        res = self.id_resolver.search_xrefs(xref, return_seq_modified=True)
        self.assertEqual(2, res[0])
        self.assertTrue(res[1], "expected that {} is a modified sequence".format(xref))

    def test_search_of_unchecked_xref(self):
        xref = "XA001g1.4"
        res = self.id_resolver.search_xrefs(xref, return_seq_modified=True)
        self.assertEqual(1, res[0])
        self.assertFalse(
            res[1], "expected that {} is a unchecked sequence".format(xref)
        )


class TaxonomyTest(unittest.TestCase):
    tax_input = None
    maxDiff = None  # show complete strings if test fails

    @classmethod
    def setUpClass(cls):
        path = find_path_to_test_db("TestDb.h5")
        h5 = tables.open_file(path)
        cls.tax_input = h5.root.Taxonomy.read()
        h5.close()

    def setUp(self):
        self.tax = Taxonomy(self.tax_input)

    def test_parents(self):
        lin = [x["Name"].decode() for x in self.tax.get_parent_taxa(284811)]
        self.assertEqual(
            lin,
            [
                "Ashbya gossypii (strain ATCC 10895 / CBS 109.51 / FGSC 9923 / NRRL Y-1056)",
                "Eremothecium",
                "Saccharomycetaceae",
                "Saccharomycetales",
                "Saccharomycetes",
                "Saccharomycotina",
                "saccharomyceta",
                "Ascomycota",
                "Dikarya",
                "Fungi",
                "Opisthokonta",
                "Eukaryota",
            ],
        )

    def test_newick(self):
        member = frozenset(
            [
                self.tax._taxon_from_numeric(x)["Name"]
                for x in self.tax.tax_table["NCBITaxonId"]
            ]
        )
        phylo = self.tax.get_induced_taxonomy(member, collapse=True)
        expected = "(((Ashbya gossypii [strain ATCC 10895 / CBS 109.51 / FGSC 9923 / NRRL Y-1056],Saccharomyces cerevisiae [strain ATCC 204508 / S288c])Saccharomycetaceae,Schizosaccharomyces pombe [strain 972 / ATCC 24843])Ascomycota,Plasmodium falciparum [isolate 3D7])Eukaryota;"
        expected = expected.replace(" ", "_")
        self.assertEqual(expected, phylo.newick())

    def test_phylogeny(self):
        member = frozenset(
            [
                self.tax._taxon_from_numeric(x)["Name"]
                for x in self.tax.tax_table["NCBITaxonId"]
            ]
        )
        phylo = self.tax.get_induced_taxonomy(member, collapse=True)
        expected = {
            "id": 2759,
            "name": "Eukaryota",
            "children": [
                {
                    "id": 4890,
                    "name": "Ascomycota",
                    "children": [
                        {
                            "id": 4893,
                            "name": "Saccharomycetaceae",
                            "children": [
                                {
                                    "id": 284811,
                                    "name": "Ashbya gossypii (strain ATCC 10895 / CBS 109.51 / FGSC 9923 / NRRL Y-1056)",
                                },
                                {
                                    "id": 559292,
                                    "name": "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)",
                                },
                            ],
                        },
                        {
                            "id": 284812,
                            "name": "Schizosaccharomyces pombe (strain 972 / ATCC 24843)",
                        },
                    ],
                },
                {"id": 36329, "name": "Plasmodium falciparum (isolate 3D7)"},
            ],
        }
        self.assertEqual(expected, phylo.as_dict())

    def test_induced_tax_simple_subtree(self):
        members = [559292, 284811]
        phylo = self.tax.get_induced_taxonomy(members)
        expected = {
            "id": 0,
            "name": "LUCA",
            "children": [
                {
                    "id": 284811,
                    "name": "Ashbya gossypii (strain ATCC 10895 / CBS 109.51 / FGSC 9923 / NRRL Y-1056)",
                },
                {
                    "id": 559292,
                    "name": "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)",
                },
            ],
        }
        self.assertEqual(expected, phylo.as_dict())

    def test_induced_tax_with_parents_subtree(self):
        members = [559292, 284811]
        phylo = self.tax.get_induced_taxonomy(members, augment_parents=True)
        expected = {
            "id": 4893,
            "name": "Saccharomycetaceae",
            "children": [
                {
                    "id": 284811,
                    "name": "Ashbya gossypii (strain ATCC 10895 / CBS 109.51 / FGSC 9923 / NRRL Y-1056)",
                },
                {
                    "id": 559292,
                    "name": "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)",
                },
            ],
        }
        self.assertEqual(expected, phylo.as_dict())

    def test_subtree_tax(self):
        clade = "Fungi"
        subtax = self.tax.get_subtaxonomy_rooted_at(clade)
        self.assertIn(4893, subtax.tax_table["NCBITaxonId"])
        self.assertNotIn(36329, subtax.tax_table["NCBITaxonId"])

    def test_taxid_of_extent_genomes(self):
        taxids = set(self.tax.get_taxid_of_extent_genomes())
        self.assertEqual({284811, 559292, 284812, 36329}, taxids)

    def test_induced_tax_with_gaps(self):
        pass


class TaxonomyTestInternalLevelSpecies(unittest.TestCase):
    """This test asserts that species A, which is also an internal level
    of the taxonomy is still reported as an external species via a
    <sciname> (disambiguate <code>) leaf."""

    taxtab = numpy.array(
        [(10, 0, b"Root"), (20, 10, b"Outgroup"), (30, 10, b"A"), (40, 30, b"B")],
        dtype=tables.dtype_from_descr(tablefmt.TaxonomyTable),
    )

    def setUp(self):
        patcher = mock.patch("pyoma.browser.models.Genome")
        self.addCleanup(patcher.stop)
        genome = patcher.start()
        type(genome).uniprot_species_code = mock.PropertyMock(return_value="HELLO")
        self.tax = Taxonomy(self.taxtab, genomes={20: genome, 30: genome, 40: genome})

    def test_newick(self):
        exp = "(Outgroup,(B,A_[disambiguate_HELLO])A)Root;"
        self.assertEqual(exp, self.tax.newick())

    def test_phyloxml(self):
        ins = (
            b"<scientific_name>A (disambiguate HELLO)</scientific_name>",
            b"<scientific_name>A</scientific_name>",
            b"<scientific_name>B</scientific_name>",
            b"<scientific_name>Outgroup</scientific_name>",
        )
        outs = (
            b"<scientific_name>B (disambiguate",
            b"<scientific_name>Outgroup (disambiguate",
        )
        res = self.tax.as_phyloxml()
        for part in ins:
            self.assertIn(part, res)
        for part in outs:
            self.assertNotIn(part, res)

    def test_dict_repr(self):
        res = self.tax.as_dict()
        self.assertDictEqual(
            {
                "name": "Root",
                "id": 10,
                "children": [
                    {"name": "Outgroup", "id": 20, "code": "HELLO"},
                    {
                        "name": "A",
                        "id": 30,
                        "children": [
                            {"name": "B", "id": 40, "code": "HELLO"},
                            {
                                "name": "A (disambiguate HELLO)",
                                "id": 30,
                                "code": "HELLO",
                            },
                        ],
                    },
                ],
            },
            res,
        )

    def test_induced_subtree_retains_internal_species(self):
        phylo = self.tax.get_induced_taxonomy([20, 40], augment_parents=True)
        self.assertIn(30, phylo.tax_table["NCBITaxonId"])


def _get_taxtab():
    # a random sample of 30 species from a production OMA release
    taxtab = numpy.array(
        [
            (2, 0, b"Bacteria"),
            (633, 1649845, b"Yersinia pseudotuberculosis"),
            (1224, 2, b"Proteobacteria"),
            (1236, 1224, b"Gammaproteobacteria"),
            (2759, 0, b"Eukaryota"),
            (5204, 451864, b"Basidiomycota"),
            (5334, 5204, b"Schizophyllum commune"),
            (5660, 2759, b"Leishmania braziliensis"),
            (5786, 2759, b"Dictyostelium purpureum"),
            (5794, 2698737, b"Apicomplexa"),
            (5807, 5794, b"Cryptosporidium parvum"),
            (5864, 5794, b"Babesia"),
            (5866, 5864, b"Babesia bigemina"),
            (6412, 33317, b"Helobdella robusta"),
            (6960, 33317, b"Hexapoda"),
            (7263, 6960, b"Drosophila arizonae"),
            (9813, 32524, b"Procavia capensis"),
            (13735, 32524, b"Pelodiscus sinensis"),
            (28211, 1224, b"Alphaproteobacteria"),
            (29760, 2759, b"Vitis vinifera"),
            (32524, 33213, b"Amniota"),
            (33154, 2759, b"Opisthokonta"),
            (33213, 33154, b"Bilateria"),
            (33317, 33213, b"Protostomia"),
            (37360, 2698737, b"Plasmodiophora brassicae"),
            (55529, 2759, b"Guillardia theta"),
            (91347, 1236, b"Enterobacterales"),
            (119060, 1224, b"Burkholderiaceae"),
            (158441, 6960, b"Folsomia candida"),
            (164328, 2698737, b"Phytophthora ramorum"),
            (214092, 1649845, b"Yersinia pestis CO92"),
            (
                264198,
                119060,
                b"Cupriavidus pinatubonensis (strain JMP 134 / LMG 1197)",
            ),
            (
                273123,
                633,
                b"Yersinia pseudotuberculosis serotype I (strain IP32953)",
            ),
            (
                284812,
                451864,
                b"Schizosaccharomyces pombe (strain 972 / ATCC 24843)",
            ),
            (
                323261,
                1236,
                b"Nitrosococcus oceani (strain ATCC 19707 / BCRC 17464 / NCIMB 11848 / C-107)",
            ),
            (338187, 1236, b"Vibrio campbellii (strain ATCC BAA-1116 / BB120)"),
            (347255, 28211, b"Rickettsia africae (strain ESF-5)"),
            (353152, 5807, b"Cryptosporidium parvum (strain Iowa II)"),
            (
                355278,
                2,
                b"Leptospira biflexa serovar Patoc (strain Patoc 1 / Ames)",
            ),
            (391774, 1224, b"Desulfovibrio vulgaris subsp. vulgaris (strain DP4)"),
            (451864, 33154, b"Dikarya"),
            (
                502801,
                633,
                b"Yersinia pseudotuberculosis serotype IB (strain PB1/+)",
            ),
            (578458, 5334, b"Schizophyllum commune (strain H4-8 / FGSC 9210)"),
            (626418, 119060, b"Burkholderia glumae (strain BGR1)"),
            (640131, 91347, b"Klebsiella variicola (strain At-22)"),
            (
                713600,
                91347,
                b"Buchnera aphidicola subsp. Acyrthosiphon pisum (strain JF98)",
            ),
            (
                754035,
                28211,
                b"Mesorhizobium australicum (strain HAMBI 3006 / LMG 24608 / WSM2073)",
            ),
            (905079, 55529, b"Guillardia theta (strain CCMP2712)"),
            (1133968, 5864, b"Babesia microti (strain RI)"),
            (1230383, 5204, b"Malassezia sympodialis (strain ATCC 42132)"),
            (1649845, 91347, b"Yersinia pseudotuberculosis complex"),
            (2698737, 2759, b"Sar"),
        ],
        dtype=tables.dtype_from_descr(tablefmt.TaxonomyTable),
    )
    return taxtab


class LucaBasedTaxonomyTests(unittest.TestCase):
    def setUp(self):
        taxtab = _get_taxtab()
        self.tax = Taxonomy(taxtab)
        self.nr_species = 30

    def test_root_is_luca(self):
        self.assertEqual(self.tax._get_root_taxon()["Name"], b"LUCA")

    def test_search_luca(self):
        node = self.tax.get_taxnode_from_name_or_taxid("LUCA")
        self.assertEqual(len(node), 1)
        self.assertEqual(node["Name"], b"LUCA")

    def test_extend_genomes(self):
        self.assertEqual(len(self.tax.get_taxid_of_extent_genomes()), self.nr_species)

    def test_search_approx_luca(self):
        self.assertIn("LUCA", [x[1] for x in self.tax.approx_search("luca")])

    def test_numeric_luca_is_found(self):
        self.assertEqual(self.tax.get_taxnode_from_name_or_taxid(0)["Name"], b"LUCA")

    def test_subtaxonomy_does_not_contain_luca_anymore(self):
        # subset of only eukaryotes
        subtax = self.tax.get_induced_taxonomy(
            [451864, 55529, 158441], collapse=True, augment_parents=True
        )
        self.assertNotIn(b"LUCA", subtax._get_root_taxon()["Name"])

    def test_as_newick_works(self):
        self.assertTrue(self.tax.newick().endswith("LUCA;"), self.tax.newick()[-30:])


class LucaWithNegTaxIDGenomeTaxonomyTests(LucaBasedTaxonomyTests):
    def setUp(self):
        taxtab = _get_taxtab()
        taxtab = numpy.append(
            taxtab,
            numpy.array([(-2, 451864, b"Some Random Genome")], dtype=taxtab.dtype),
        )
        self.tax = Taxonomy(taxtab)
        self.nr_species = 31

    def test_neg_taxid_species_in_newick(self):
        self.assertIn("Some_Random_Genome", self.tax.newick())


class DBMock(object):
    def __init__(self, h5):
        self.h5 = h5

    def get_hdf5_handle(self):
        return self.h5


class MockDBTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = find_path_to_test_db("TestDb.h5")
        cls.db = DBMock(tables.open_file(path))

    @classmethod
    def tearDownClass(cls):
        cls.db.get_hdf5_handle().close()


class GenomeIdResolverTest(MockDBTestCase):
    def setUp(self):
        self.OmaIdMapper = OmaIdMapper(self.db)

    def test_resolve_uniprot_code(self):
        query = "YEAST"
        g = self.OmaIdMapper.identify_genome(query)
        self.assertEqual(g["UniProtSpeciesCode"], query.encode("ascii"))

    def test_resolve_non_existing_uniprot_code(self):
        query = "XED22"
        with self.assertRaises(UnknownSpecies):
            self.OmaIdMapper.identify_genome(query)

    def test_resolve_taxon_id(self):
        expected = b"PLAF7"
        for query in (36329, b"36329", "36329"):
            g = self.OmaIdMapper.identify_genome(query)
            self.assertEqual(
                g["UniProtSpeciesCode"],
                expected,
                "failed for {} (type {})".format(query, type(query)),
            )

    def test_resolve_nonexisting_code(self):
        with self.assertRaises(UnknownSpecies):
            self.OmaIdMapper.identify_genome(2)

    def test_approx_search_genome(self):
        query = "sacero cervesa"
        expect = "YEAST"
        cands = self.OmaIdMapper.approx_search_genomes(query)
        self.assertIn(expect, [g.uniprot_species_code for g in cands])


class TestPerGenomeMetaData(MockDBTestCase):
    def setUp(self) -> None:
        self.pg = PerGenomeMetaData(self.db.get_hdf5_handle(), "YEAST")

    def test_in_oma_groups_matches(self):
        pass


class FastMapperTester(unittest.TestCase):
    def setUp(self) -> None:
        self.fast_mapper = FastMapper(Database(find_path_to_test_db()))

    def tearDown(self) -> None:
        self.fast_mapper.db.close()

    def test_search_small_sequence(self):
        test_seq = io.StringIO(">test\nAS")
        seq_iter = SeqIO.parse(test_seq, format="fasta")
        with self.assertLogs("pyoma", level=logging.INFO) as cm:
            anno = list(self.fast_mapper.iter_projected_goannotations(seq_iter))
        self.assertIn("Skipping", "\n".join(cm.output))
        self.assertEqual(0, len(anno))

    def check_mapped_seqs(self, seqs, way="Exact"):
        with self.assertLogs("pyoma", level=logging.INFO) as cm:
            anno = list(self.fast_mapper.iter_projected_goannotations(seqs))
        for a in anno:
            with_ = self.fast_mapper.db.id_mapper["OMA"].map_entry_nr(
                int(a["DB_Object_ID"])
            )
            self.assertEqual(f"{way}:{with_}".split(":"), a["With"].split(":")[:2])
        self.assertGreaterEqual(len(a), 1)

    def test_search_existing_sequences(self):
        seqs = [
            SeqRecord(
                id=f"{enr}", seq=Seq(self.fast_mapper.db.get_sequence(enr).decode())
            )
            for enr in range(1, 5)
        ]
        self.check_mapped_seqs(seqs, way="Exact")

    def test_search_inexact_search_seq(self):
        seqs = []
        for en in range(1, 5):
            s = self.fast_mapper.db.get_sequence(en).decode()
            modif_seq = s.replace("A", "R")
            seqs.append(SeqRecord(id=f"{en}", seq=Seq(modif_seq)))
        self.check_mapped_seqs(seqs, way="Approx")
