import json
import os
import shutil
import tempfile
import unittest
from hashlib import md5

from pyoma.browser import tablefmt
from pyoma.browser.build.builder import DBBuilder, OmaGroupsProvider


def _write_genome_json(path, nr_entries, *, seq_len=12, offset=0):
    """Create a minimal per-genome protein-data JSON as consumed by add_proteins."""
    data = {
        "seqs": [f"{'M' * (seq_len - 1)}{i}" for i in range(nr_entries)],
        "cdna": [f"{'ATG' * seq_len}{i}" for i in range(nr_entries)],
        "de": [f"description {offset + i + 1}" for i in range(nr_entries)],
        "chrs": ["chr1" for _ in range(nr_entries)],
        "acs": [[f"AC{offset + i + 1:04d}"] for i in range(nr_entries)],
        "locs": [[[1, 3 * seq_len, 1]] for _ in range(nr_entries)],
    }
    with open(path, "w") as fh:
        json.dump(data, fh)
    return data


class OmaGroupsProviderTest(unittest.TestCase):
    """Unit tests for the genome-code / genome-name translation added to OmaGroupsProvider."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.groups_path = os.path.join(self.tmpdir, "oma_groups.json")
        with open(self.groups_path, "w") as fh:
            json.dump({"Saccharomyces cerevisiae S288C": {"1": 10, "2": 0}, "Homo sapiens": {"1": 20}}, fh)
        self.name2code = {
            "Saccharomyces cerevisiae S288C": "YEAST",
            "Homo sapiens": "HUMAN",
        }

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_lookup_translates_code_to_name(self):
        provider = OmaGroupsProvider(self.groups_path, self.name2code)
        self.assertEqual(provider.get_oma_group("YEAST", 1), 10)
        self.assertEqual(provider.get_oma_group("YEAST", 2), 0)
        self.assertEqual(provider.get_oma_group("HUMAN", 1), 20)

    def test_unknown_genome_code_returns_zero(self):
        provider = OmaGroupsProvider(self.groups_path, self.name2code)
        self.assertEqual(provider.get_oma_group("FOOBA", 1), 0)

    def test_missing_entry_nr_returns_zero(self):
        provider = OmaGroupsProvider(self.groups_path, self.name2code)
        self.assertEqual(provider.get_oma_group("YEAST", 99), 0)

    def test_without_name2code_falls_back_to_raw_key(self):
        # backward compatible behaviour: data keyed directly by whatever
        # identifier is passed to get_oma_group (e.g. by genome code)
        with open(self.groups_path, "w") as fh:
            json.dump({"YEAST": {"1": 10}}, fh)
        provider = OmaGroupsProvider(self.groups_path)
        self.assertEqual(provider.get_oma_group("YEAST", 1), 10)


class AddProteinsTest(unittest.TestCase):
    """Integration test for DBBuilder.add_proteins, in particular that the OmaGroup
    is correctly resolved via OmaGroupsProvider when genome codes differ from the
    genome names used as keys in the oma-groups json file."""

    GENOMES = [
        # code, sci_name, taxid, tot_entries, entry_off
        ("YEAST", "Saccharomyces cerevisiae S288C", 4932, 2, 0),
        ("HUMAN", "Homo sapiens", 9606, 1, 2),
    ]

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.db_path = os.path.join(self.tmpdir, "test.h5")
        self.name2code = {name: code for code, name, _, _, _ in self.GENOMES}

        self.code_to_file = {}
        self.genome_data = {}
        for code, _, _, tot_entries, entry_off in self.GENOMES:
            fname = os.path.join(self.tmpdir, f"{code}.json")
            self.genome_data[code] = _write_genome_json(fname, tot_entries, offset=entry_off)
            self.code_to_file[code] = fname

        oma_groups_path = os.path.join(self.tmpdir, "oma_groups.json")
        with open(oma_groups_path, "w") as fh:
            json.dump(
                {
                    "Saccharomyces cerevisiae S288C": {"1": 42, "2": 0},
                    "Homo sapiens": {"1": 7},
                },
                fh,
            )
        self.oma_group_provider = OmaGroupsProvider(oma_groups_path, self.name2code)
        self.expected_oma_group = {1: 42, 2: 0, 3: 7}

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _create_genome_table(self, db):
        gs_tab = db.h5.create_table("/", "Genome", tablefmt.GenomeTable, expectedrows=len(self.GENOMES))
        for code, sci_name, taxid, tot_entries, entry_off in self.GENOMES:
            row = gs_tab.row
            row["NCBITaxonId"] = taxid
            row["UniProtSpeciesCode"] = code
            row["TotEntries"] = tot_entries
            row["TotAA"] = tot_entries * 20
            row["EntryOff"] = entry_off
            row["SciName"] = sci_name
            row["IsPolyploid"] = False
            row.append()
        gs_tab.flush()
        return gs_tab

    def test_add_proteins_assigns_oma_group_via_genome_code_translation(self):
        with DBBuilder(self.db_path, mode="write") as db:
            self._create_genome_table(db)
            db.add_proteins(self.code_to_file, self.oma_group_provider, xref_collector=None)

            entry_tab = db.h5.get_node("/Protein/Entries")
            self.assertEqual(len(entry_tab), 3)

            seen = {}
            for row in entry_tab:
                seen[int(row["EntryNr"])] = int(row["OmaGroup"])
            self.assertEqual(seen, self.expected_oma_group)

    def test_add_proteins_writes_matching_sequence_hashes(self):
        with DBBuilder(self.db_path, mode="write") as db:
            self._create_genome_table(db)
            db.add_proteins(self.code_to_file, self.oma_group_provider, xref_collector=None)

            entry_tab = db.h5.get_node("/Protein/Entries")
            seq_buf = db.h5.get_node("/Protein/SequenceBuffer")
            for row in entry_tab:
                seq = seq_buf[row["SeqBufferOffset"] : row["SeqBufferOffset"] + row["SeqBufferLength"] - 1].tobytes()
                self.assertEqual(md5(seq).hexdigest(), row["MD5ProteinHash"].decode())


if __name__ == "__main__":
    unittest.main()
