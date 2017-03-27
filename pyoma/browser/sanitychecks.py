""" Sanity check for each release of OMA hdf5 databases

Description: Natasha leads

.. moduleauthor:: Natasha Glover <gfhtk>, Henning Redestig <gbfjc>, Jiao Long <ggqcq>

"""

from collections import Counter
from functools import reduce
import re
import os
import json

import tables
from scipy import stats
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
import numpy as np
import ggplot
from prettytable import PrettyTable as PT
from BCBio import GFF
import matplotlib.pyplot as plt
from Bio import SeqIO

from bcsoma import omadatautils as od

darwin_browser_share = '/tools/bioinfo/app/OMA_pipeline-{0}/data/Browser/'

darwin_browserdata_path= '/tools/bioinfo/app/OMA_pipeline-{0}/data/Browser/{1}/data/OmaServer.h5'

class SanitySession(object):
    """A class for performing descriptive analysis of a OMA browser database"""
    
    def __init__(self, oma_instance = None, release = None):
        """Constructor
        :param oma_ins: OMA instance i.e. production/qa/dev
        :param realease: OMA OmaServer.h5 db release e.g. All.May2016
        """
        self.ins = oma_instance
        self.release = release
        self.db_path = None
        if os.path.exists(darwin_browser_share.format(self.ins,self.release)):
            self.db_path = darwin_browserdata_path.format(self.ins,self.release)
        else:
            raise IOError(2,'OmaServer.h5 db does not exist')
        self.db = self.read_oma_db()
        self.species = self.get_species().keys()
        self.entries_table = self.db.root.Protein.Entries
        self.hoglevel_table = self.db.root.HogLevel
        self.omagroups = self.get_omagroups()
        self.all_hogs = self.get_all_hogs()
        self.all_hog_lvls = self.get_all_hog_lvls()
        self.all_lvls = self.get_all_lvls() 
        #self.describe = stats_omagroup()
        #self.omagroups_df = genes_per_omagroup()

    def read_oma_db(self):
        return tables.open_file(self.db_path, mode = "r")

    def close_oma_db(self):
        self.db.close()

    def get_species(self):
        """returns a dictionary genomes that the key is OMA 5 letters code and the value
        is taxon identifier of each species in this release

        :param taxonid: returns OMA five letter code if True, otherwise returns tax id 
        """
        summaries_file = os.path.join(os.path.dirname(self.db_path),'Summaries.drw')
        genomes  = {}
        with open(summaries_file) as summaries:
            for x in summaries:
                if re.match(r'^GenomeSummaries\[.+', x):
                    five_letter = re.findall(r'([A-Z0-9]{5})', x)[0]
                if re.match(r'^GenomeSummary\(.+', x):
                    taxid = re.findall(r'<TAXONID>([0-9]+)</TAXONID>', x)[0]
                    genomes[five_letter] = taxid
        return genomes
                            
    def get_omagroups(self):
        """get the number of genes in each oma group and returns a dictionary/Counter
        object that key: oma group id, value: number of genes in this oma group 
        """
        omagroups = Counter([x['OmaGroup'] for x in self.entries_table.where('OmaGroup!=0')])
        return omagroups

    def get_all_lvls(self):
        """get all taxonomic levels"""
        return [x['Level'] for x in self.hoglevel_table.iterrows()]
    
    def get_all_hogs(self):
        """get the number of genes in HOGs at all taxonomic clades and returns a dictionary/Counter
        object that key: HOG id, value: number of genes in this HOG at all taxonomic clades 
        """
        condition = ''
        hogs = Counter([x['OmaHOG'] for x in
                self.entries_table.where('OmaHOG!=condition')])
        return hogs

    def get_all_hog_lvls(self):
        """get the number of genes in a HOG at each taxonomic clades and returns a dictionary/Counter
        object that key: taxonomic clade, value: number of genes in a HOG at this taxonomic clades 
        """
        levels = Counter([x['Level'] for x in self.hoglevel_table.iterrows()])
        return levels 
    
#    def stats_omagroup(self):
#        """calculates descriptive statistics about the OMA groups"""
#        
#        val = list(self.omagroups.values())
#        describe = stats.describe(val)
#        result = {}
#        result['release'] = str(self.release)
#        result['Nb genes in OmaGroups'] = sum(val)
#        result['Nb genes NOT in OmaGroups'] = self.entries_table.nrows - sum(val)
#        result['Nb OmaGroups'] = describe[0]
#        result['min, max nb genes per OmaGroup'] = describe[1]
#        result['mean nb genes per OmaGroup'] = describe[2]
#        return result

#    def stats_all_hog(self):
#        """Calculates descriptive statistics on all HOGs"""
#        val = list(self.all_hogs.values())
#        describe = stats.describe(val)
#        result = {}
#        result['release'] = self.release
#        result['Nb HOGs'] = describe[0]
#        result['Nb genes in HOGs'] = sum(val)
#        result['% of total genes in HOGs'] = float(sum(val) /
#                                                   self.entries_table.nrows)
#        result['mean nb genes per HOG'] = describe[2]
#       result['min, max nb genes per HOG'] = describe[1]
#       return result

#    def stats_hoglevels(self):
#        "calculates the count distribution for each taxonomic level"
#        val = list(self.all_hog_lvls.values())
#        describe = stats.describe(val)
#        result = {}
#        result['Nb HOGs for each level'] = self.all_hog_lvls
#        result['Nb levels in total'] = describe[0]
#        result['Level with most HOGs'] = max(self.all_hog_lvls, key=self.all_hog_lvls.get)
#        result['Level with least HOGs'] = min(self.all_hog_lvls, key=self.all_hog_lvls.get)
#        result['mean nb HOGs per level'] = describe[2]
#        result['release'] = str(self.release)
#        return result                                               
    
    def genes_per_omagroup(self):
        """get the numbers of gene per oma group as a data frame"""
        result = pd.DataFrame(pd.Series(self.omagroups), columns=['count'])
        result.index.names = ['omagroup']
        result['release'] = self.release
        return result

    def genes_per_hog(self):
        """get the numbers of gene per hog as a data frame"""
        result = pd.DataFrame(pd.Series(self.all_hogs), columns=['count'])
        result.index.names = ['hog']
        result['release'] = self.release
        return result

    def R_genes_per_omagroup(self):
        """This function prints the number of genes for each OMA group
        as a csv formatted long string...
        """
        result = ""
        for k in self.omagroups:
            v = self.omagroups[k]
            result += self.release + "\t" + str(k) + "\t" + str(v) + "\n"
        return (result)

    def R_genes_per_hog(self):
        """prints the number of genes for each HOG"""
        result = ""
        for k in self.all_hogs:
            v = self.all_hogs[k]            
            result += self.release + "\t" + str(k.decode("utf-8")) + "\t" + str(v) + "\n"
        return result

    def num_hogs_per_level(self):
        """get the number of hogs for a list of taxa levels"""
        result = ""
        for taxa in self.all_lvls:
            result += str(taxa.decode("utf-8"))+"\t"+ self.release+"\t"+str(self.all_hog_lvls.get(taxa))+"\n"
            result = re.sub(r'None', "", result)
        return result

def get_allopolypoid_genomes():
    """get allopolypoid genomes in OMA"""
    allopolypoid_genomes = ['WHEAT', 'GOSHI', 'BRANA']
    return allopolypoid_genomes
    
def homoeolog_checks(**kwargs):
    columns = ['species','release','Nb homoeolog pairs']
    rel1 = rel2 = None
    try:
        rel1 = kwargs.get('rel1')
        rel2 = kwargs.get('rel2')
    except ValueError:
        print ('homeolog_checks() expects two releases')
    if (not valid_release(rel1.ins,rel1.release)) or (not valid_release(rel2.ins,rel2.release)):
        raise OSError(2, 'Release does not exist')   
    species_list = None
    if kwargs.get('species'):
        species_list = kwargs.get('species')
        if (species_list not in rel1.species) or (species_list not in rel2.species):
            raise ValueError('The given two OMA releases do not cover all the query species')
    else:
        species_list = get_allopolypoid_genomes()
    print(species_list)
    homo_df = pd.DataFrame(columns=columns)
    for spec in species_list:
        for rel in (rel1,rel2):
            child = rel.db.root.PairwiseRelation._f_get_child(spec)
            homoeologous_pairs = child._f_list_nodes()[1].read_where('RelType == 5')
            hp_row_ins = pd.DataFrame({'species': [spec],
                                       'release': [rel.release],
                                       'Nb homoeolog pairs': [len(homoeologous_pairs)/2]})
            homo_df = homo_df.append(hp_row_ins,ignore_index = True)
    fig = plt.figure()
    ax = sns.pointplot(x="release", y="Nb homoeolog pairs", hue='species',data=homo_df)
    #plt.legend(loc='upper left')
    fig.add_axes(ax)
    fig.savefig('homoeologs_check.png')

def get_sorted_rels(oma_insance):
    """get releases from a OMA instance sorted by the time of last modification
    :param oma_ins: Bayer inhouse OMA instance e.g. 'production', 'qa', or 'dev'
    """
    dir_path = darwin_browser_share.format(oma_insance)
    releases = [s for s in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path,s))]
    releases.sort(key=lambda s: os.path.getmtime(os.path.join(dir_path, s)), reverse = True)
    return releases

def valid_release(oma_instance,rel):
    """validate whether the given release name exists
    :param rel: OMA release
    """
    return os.path.exists(os.path.join(darwin_browser_share.format(oma_instance),rel))

def consistency_check(aa_file, nt_file, gff_file, gene_id_re, lvl, output):
    """consistency check on the records among amino acid sequence file, nucleotide sequence file, and 
    GFF file. Output to stdout or to a report file. 
    """
    aa_records = get_gene_list(aa_file,gene_id_re)
    nt_records = get_gene_list(nt_file,gene_id_re)
    gff_records = get_protein_coding_genes(get_gff_records(gff_file), gene_id_re, lvl)
    intersection_records = reduce(set.intersection,[aa_records,nt_records,gff_records])
    with open(output,'w') as out:
        out.write('The number of records in amino acid sequence file: {}\n'.format(len(aa_records)))
        out.write('The number of records in nucleotide sequence file: {}\n'.format(len(nt_records)))
        out.write('The number of records in GFF file: {}\n'.format(len(gff_records)))
        out.write('ERROR: Pay attention to the following records\n')
        for rec in intersection_records:
            out.write('{}\n'.format(rec))
        
    assert(aa_records == nt_records == gff_records),'records are not matched among genome resource files'
    return gff_records
    
def get_gene_id(record_header, gene_id_re=''):
    """get only gene identifier by trimming off sub string after first column, 
    spaces, and mRNA tail from the header line of a records in a fasta file.
    :param gene_id_re: trim off pattern match   
    """
    return(re.sub(gene_id_re,'',record_header))

def get_gene_list(seq_file, gene_id_re):
    """retrieve all records in a fasta file
    :param seq_file: sequence file e.g. aa/nt file
    """
    return([get_gene_id(record.id,gene_id_re) for record in SeqIO.parse(seq_file,"fasta")])

def get_gff_records(gff_file):
    """stores all gene records in a gff file into a list
    :param gff_file: gff format file
    """
    gff_records = []
    with open(gff_file) as gff:
        for record in GFF.parse(gff):
            gff_records.append(record)
    return gff_records
def get_protein_coding_genes(gff_records,gene_id_re,lvl):
    """retrives only protein-coding genes identifiers from a GFF file
    :gff_records: gene records list 
    """
    prot_coding_gene = []
    for chromosome in gff_records:
        for gene in chromosome.features:
            if gene.type == lvl:
                gene_id = re.sub(gene_id_re, '',  gene.id)
                if gene_id not in prot_coding_gene:
                    prot_coding_gene.append(gene_id)
    return prot_coding_gene

class OMASpeciesDB(object):

    def __init__(self):
        """Constructor"""
        self.db = os.path.join(os.environ['DARWIN_GENOMES_PATH'],'inhouse_specs.json')
        self.specs_data = self.get_species_info()
        self.species = self.get_species()
        
    def get_species_info(self):
        """get species records"""
        with open(self.db) as spec_db:
            specs_info = json.load(spec_db)
        return specs_info

    def get_species(self):
        return([ spec_rec['species'] for spec_rec in self.specs_data])
    
    def add_species(self, spec, aa, nt, gff, nr_recs):
        """add a species record (Be aware of duplicate records)
        :param spec: five code represents a species to add
        :param aa: path to protein sequence file
        :param nt: path to nucleotide sequence file
        :param gff: path to GFF file
        :nr_recs: the number of protein records in this release
        """
        if spec_ not in self.species:
            with open(self.db, 'w') as spec_db:
                specs_ins = {}
                specs_ins['species'] = spec_
                specs_ins['aa'] = aa
                specs_ins['nt'] = nt
                specs_ins['gff'] = gff
                specs_ins['nr_recs'] = nr_recs
                self.specs_data.append(specs_ins)
                json.dump(self.specs_data, spec_db, sort_keys = True, indent=4)
        else:
            raise ValueError('Genome resource of {} is available in the current db'\
                             .format(spec))
    
        
if __name__ == "__main__":
    ####################################################################
    #oma_ins = 'dev'
    #rels = get_sorted_rels(oma_ins)
    #cur_rel = rels[0]
    #pre_rel = rels[1] 
    #cur_rel_ss = SanitySession(oma_instance = oma_ins, release = cur_rel)
    #pre_rel_ss = SanitySession(oma_instance = oma_ins, release = pre_rel)
    #homoeolog_checks(rel1 = cur_rel_ss, rel2 = pre_rel_ss)
    ####################################################################
    newDB = OMASpeciesDB()
    path = '/tools/bioinfo/app/OMA_pipeline-20140123/data/OMA/genomes_incoming_inhouse/'
    #Genome resources of CAPAN(pepper) and WHEAT integrated into OMA are from public
    #CAPAN public resource (checked)
    #BRAJU genome assembly is littered with transposons (need to filter out transposons)
    #BRANA adopted from PLAZA build Darmor jira: GUS-29869 (needs to fix)
    #CUCSA 

    inhouse_specs = [ {'five_letter':'BRANA', 'geneid_re':'', 'gene_lvl' 'mRNA',\
                      {'five_letter':'GOSRA', 'geneid_re':},\
                      {'five_letter':'BRAOA', 'trim_dna_id_re':'^[^:]+:', 'trim_protein_id_re':'^[^:]+:', 'geneid_re':'\.\d+$'},
                      {'five_letter':'TRITU', 'trim_dna_id_re':['^[^:]+:'], 'trim_protein_id_re':['^[^:]+:', '/.+$'], 'trim_protein_desc_re':[".*"]},\
                      {'five_letter':'CUCSA', 'geneid_re':'\.\d$'},\
                      {'five_letter':'BRARP', 'trim_protein_id_re':['-P$'], 'geneid_re':'\.\d+$', 'gff_mrna_name':'ID','gff_mrna_id_re':'^transcript:'},\
                      {'five_letter':'CITLA', 'gff_mrna_name':'ID'},\
                      {'five_letter':'GOSHI', 'trim_protein_id_re':['-P$'], 'geneid_re':'^evm..U.', 'gff_mrna_name':'ID'},\
                      {'five_letter':'CUCME', 'geneid_re':'T[0-9]+$', 'gff_mrna_name':'ID', 'mrna_feature':'transcript'}]
    
    for spec in inhouse_species:
        spec_aa_file, spec_nt_file, spec_gff_file = map(lambda x: os.path.join(path,spec,'data'),\
                                        ['{}.aa.fasta'.format(spec['five_letter']),\
                                         '{}.nt.fasta'.format(spec['five_letter']),\
                                         '{}.gff3'.format(spec['five_letter'])])
        gene_id_re = spec['geneid_re']
        lvl = spec['gff_mrna_id_re'] #top level normally either gene or mRNA
        nr_records = consistency_check(spec_aa_file, spec_nt_file, spec_gff_file, gene_id_re, lvl, output)
        newDB.add_species(spec,spec_aa_file,spec_nt_file,spec_gff_file, nr_records)
