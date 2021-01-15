import tables
import pandas as pd
import os
import re
from tqdm import tqdm
from collections import Counter

oma_browserdata_path = "{0}/{1}/data/OmaServer.h5"


class SanitySession(object):
    """A class for performing descriptive analysis of a OMA browser database"""

    def __init__(self, oma_browser_dir=None, release=None):
        """
        :param oma_ins: OMA instance i.e. production/qa/dev
        :param realease: OMA OmaServer.h5 db release e.g. All.May2016
        """
        self.release = release
        self.db_path = None
        if os.path.exists(oma_browserdata_path.format(oma_browser_dir, self.release)):
            self.db_path = oma_browserdata_path.format(oma_browser_dir, self.release)
        else:
            raise IOError(2, "OmaServer.h5 db does not exist")

        self.h5_handle = self.read_oma_db()
        self.species = self.get_species2()

        self.entries_table = self.h5_handle.root.Protein.Entries
        self.genome_table = self.h5_handle.root.Genome
        self.hoglevel_table = self.h5_handle.root.HogLevel
        self.omagroups = self.get_omagroups()
        self.all_hogs = self.get_all_hogs()
        self.all_hog_lvls = self.get_all_hog_lvls()
        self.all_lvls = self.get_all_lvls()
        self.xref_df = self.get_number_of_xrefs_per_species()
        print(self.release, "done")

    def read_oma_db(self):
        return tables.open_file(self.db_path, mode="r")

    def close_oma_db(self):
        self.h5_handle.close()

    def get_species2(self):
        """Uses the h5 file to get species"""
        genome_df = pd.DataFrame(self.h5_handle.root.Genome.read())
        genomes = genome_df["UniProtSpeciesCode"]
        genomes = [x.decode("utf-8") for x in genomes]
        return genomes

    def get_omagroups(self):
        """get the number of genes in each oma group and returns a dictionary/Counter
        object that key: oma group id, value: number of genes in this oma group
        """
        omagroups = Counter(
            (x["OmaGroup"] for x in self.entries_table.where("OmaGroup!=0"))
        )
        return omagroups

    def get_all_lvls(self):
        """get all taxonomic levels"""
        return [x["Level"] for x in self.hoglevel_table.iterrows()]

    def get_all_hogs(self):
        """get the number of genes in HOGs at all taxonomic clades and returns a dictionary/Counter
        object that key: HOG id, value: number of genes in this HOG at all taxonomic clades
        """
        hogs = Counter()
        dot_pat = re.compile(b"\.")
        for prot in self.entries_table.where("OmaHOG != b''"):
            hogid = prot["OmaHOG"]
            hogs[hogid] += 1
            for m in dot_pat.finditer(hogid):
                hogs[hogid[: m.start()]] += 1
        return hogs

    def get_all_hog_lvls(self):
        """get the number of genes in a HOG at each taxonomic clades and returns a dictionary/Counter
        object that key: taxonomic clade, value: number of genes in a HOG at this taxonomic clades
        """
        levels = Counter([x["Level"] for x in self.hoglevel_table.iterrows()])
        return levels

    def genes_per_omagroup(self):
        """get the numbers of gene per oma group as a data frame"""
        result = pd.DataFrame(pd.Series(self.omagroups), columns=["count"])
        result.index.names = ["omagroup"]
        result["release"] = self.release
        return result

    def genes_per_hog(self):
        """get the numbers of gene per hog as a data frame"""
        result = pd.DataFrame(pd.Series(self.all_hogs), columns=["count"])
        result.index.names = ["hog"]
        result["release"] = self.release
        return result

    def num_hogs_per_level(self):
        """get the number of hogs for a list of taxa levels"""
        result = ""
        for taxa in self.all_lvls:
            result += (
                str(taxa.decode("utf-8"))
                + "\t"
                + self.release
                + "\t"
                + str(self.all_hog_lvls.get(taxa))
                + "\n"
            )
            result = re.sub(r"None", "", result)
        return result

    def get_number_of_xrefs_per_species(self):
        """load all the xrefs we have of each type and species"""
        xreftab = self.h5_handle.root.XRef
        enum = xreftab.get_enum("XRefSource")
        genomes = self.h5_handle.root.Genome.read()
        sorter = genomes.argsort(order=("EntryOff"))

        def map_to_species(entry_nrs):
            idx = genomes["EntryOff"].searchsorted(
                entry_nrs - 1, side="right", sorter=sorter
            )
            return genomes[idx - 1]["UniProtSpeciesCode"]

        cnts = Counter()
        chunksize = xreftab.chunkshape[0]
        for start in tqdm(
            range(0, len(xreftab), chunksize), desc="processing xref-chunks"
        ):
            chunk = xreftab.read(start, start + chunksize)
            orgs = map_to_species(chunk["EntryNr"])
            for org, xref_instance in zip(orgs, chunk):
                cnts[(org, enum(xref_instance["XRefSource"]))] += 1
        df = pd.DataFrame(
            [(org.decode(), source, counts) for (org, source), counts in cnts.items()],
            columns=["Species", "Source", "Counts"],
        )
        return df
