#!/usr/bin/env python3

import csv
import itertools
import operator
import pandas
import numpy
import pyoma.browser.db as pydb
import pyoma.browser.models as models
import numpy
import gzip
import os
import logging

logger = logging.getLogger('genevestigator-export')


class ColumnRenameDictWriter(csv.DictWriter):
    def __init__(self, fh, fieldnames, fieldlabels, **kwargs):
        super(ColumnRenameDictWriter, self).__init__(fh, fieldnames, **kwargs)
        self.fieldlabels = fieldlabels

    def writeheader(self):
        header = dict(zip(self.fieldnames, self.fieldlabels))
        self.writerow(header)


class GenomeFilter(object):
    def __init__(self, db, keep):
        genomes = set([])
        for node in keep:
            t = db.tax.get_subtaxonomy_rooted_at(node)
            extent_ncbi = t.get_taxid_of_extent_genomes()
            extent_genomes = [models.Genome(db, db.id_mapper['OMA'].genome_from_taxid(z)) for z in extent_ncbi]
            logger.info("{} maps to {}".format(node, [g.uniprot_species_code for g in extent_genomes]))
            genomes = genomes.union(set(extent_ncbi))
        self.genomes = sorted([models.Genome(db, db.id_mapper['OMA'].genome_from_taxid(z)) for z in genomes],
                              key=operator.attrgetter('entry_nr_offset'))
        self.offsets = numpy.array([(x.entry_nr_offset, x.nr_entries+x.entry_nr_offset) for x in self.genomes],
                                   dtype=[('begin', 'i4'), ('end', 'i4')])
        logger.info("Filter with {} species created: {}".format(len(self.genomes), self.genomes))

    def filter_entry_nr(self, entry_nr):
        """filter a numpy array of entry_nrs or a single entry_nr whether
        they belong to a species of interest.

        Returns:
            numpy array of same size as input entry_nr with boolean values
        """
        idx = numpy.searchsorted(self.offsets['begin'], entry_nr-1, side="right")
        return (self.offsets[idx-1]['begin'] < entry_nr) & (entry_nr <= self.offsets[idx-1]['end'])

    def filter_entry_pair(self, entry_nr1, entry_nr2):
        return self.filter_entry_nr(entry_nr1) & self.filter_entry_nr(entry_nr2)


class PairsExtractor(object):
    def __init__(self, db, filter):
        self.filter = filter
        self.db = db
        self._xrefs_cache = {}

    def process(self, outfh):
        cols = ['Protein1', 'ProteinAC1', 'GeneID1', 'Protein2', 'ProteinAC2', 'GeneID2',
                'OrthologyType', 'Score', 'Distance']
        labs = (*cols[:7], 'AlignmentScore', 'PamDistance')
        csv_writer = ColumnRenameDictWriter(outfh, cols, labs,
                                    extrasaction='ignore', delimiter='\t')
        csv_writer.writeheader()
        for g1, g2 in itertools.combinations(self.filter.genomes, 2):
            for rel in self.extract_pairwise_relations(g1, g2):
                csv_writer.writerow(rel)

    def _load_xrefs(self, genome, source):
        xref_tab = self.db.db.get_node('/XRef')
        xrefs = pandas.DataFrame(xref_tab.read_where(
            '(EntryNr > {}) & (EntryNr <= {}) & (XRefSource == {})'.format(
                genome.entry_nr_offset, genome.entry_nr_offset + genome.nr_entries,
                xref_tab.get_enum('XRefSource')[source])))
        xrefs.drop(columns=['XRefSource', 'Verification'], inplace=True)
        xrefs['XRefId'] = xrefs['XRefId'].str.decode("utf-8")
        xrefs_one_ref = xrefs.drop_duplicates(subset='EntryNr', keep="first")
        return xrefs_one_ref

    def get_or_load_xrefs_df(self, g, source):
        if not (g.uniprot_species_code, source) in self._xrefs_cache:
            self._xrefs_cache[(g.uniprot_species_code, source)] = self._load_xrefs(g, source)
        return self._xrefs_cache[(g.uniprot_species_code, source)]

    def augment_with_ids(self, df, g1, g2, source, label):
        def _augment_id_of_genome(nr):
            lab = "{}{}".format(label, nr)
            on = "EntryNr{}".format(nr)
            g = g1 if nr==1 else g2
            xref = self.get_or_load_xrefs_df(g, source)
            new_df = df.merge(xref, how='left', left_on=on, right_on="EntryNr")
            new_df.drop(columns=['EntryNr'], inplace=True)
            new_df.rename(columns={'XRefId': lab}, inplace=True)
            return new_df
        df = _augment_id_of_genome(1)
        df = _augment_id_of_genome(2)
        return df

    def extract_pairwise_relations(self, g1, g2):
        tab = self.db.db.get_node('/PairwiseRelation/{}/VPairs'.format(g1.uniprot_species_code))
        df = pandas.DataFrame(tab.read_where(
            '(EntryNr2 > {}) & (EntryNr2 <= {})'.format(
                g2.entry_nr_offset, g2.entry_nr_offset + g2.nr_entries)))
        df['OrthologyType'] = df['RelType'].apply(tab.get_enum('RelType'))
        df['Protein1'] = df['EntryNr1'].apply(
            lambda enr: "{}{:05d}".format(g1.uniprot_species_code, enr - g1.entry_nr_offset))
        df['Protein2'] = df['EntryNr2'].apply(
            lambda enr: "{}{:05d}".format(g2.uniprot_species_code, enr - g2.entry_nr_offset))
        df = self.augment_with_ids(df, g1, g2, 'SourceAC', 'ProteinAC')
        df = self.augment_with_ids(df, g1, g2, 'SourceID', 'GeneID')
        yield from df.to_dict('records')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Export Pairwise Orthologs in TSV format suitable for GeneVestigator")
    parser.add_argument('--out', '-o', default='-', type=argparse.FileType('w'),
                        help='Path to output file. Defaults to stdout')
    parser.add_argument('db', help="Path to hdf5 database")
    parser.add_argument('taxon', nargs="+", help="Taxon names of species/clades to be included. "
                                                 "Can be scientific name or ncbi taxon id.")
    parser.add_argument('-v', default=0, action="count", help="increase level of verbosity to INFO / DEBUG")
    conf = parser.parse_args()
    logging.basicConfig(level=30 - 10 * min(conf.v, 2))

    db = pydb.Database(conf.db)
    genome_filter = GenomeFilter(db, conf.taxon)
    extractor = PairsExtractor(db, genome_filter)
    extractor.process(conf.out)
    conf.out.close()
