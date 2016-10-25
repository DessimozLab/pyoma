__author__ = 'admin'

import familyanalyzer as fa
import lxml.etree as etree
import os

def indent(elem, level=0):
    """
    re structure the xml tree in human readable format (pre processing before writing the tree in a file)
    :param elem:
    :param level:
    :return:
    """
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

class OrthoXMLSplitter(object):

    def __init__(self, xml_file, list_hog = None):
        self.xml_file = xml_file
        self.list_hog = list_hog
        self.cache_dir = None
        self.Etree_XML = etree.parse(self.xml_file)
        self.Etree_root = self.Etree_XML.getroot()
        self.Etree_OGs = fa.OrthoXMLQuery.getToplevelOrthologGroups(self.Etree_root)
        self.Etree_header_genes = fa.OrthoXMLQuery.getInputGenes(self.Etree_root)

    def __call__(self, cache_dir):
        self.cache_dir = cache_dir
        for og in self.Etree_OGs:
            if self.list_hog is None:
                self.write_OG(og)
            else:
                if int(og.get("id")) in self.list_hog:
                    self.write_OG(og)

    def get_generef_OG(self, Etree_OG):
        generef_els = []
        for gene_etree in Etree_OG.getiterator():
            if gene_etree.tag == "{http://orthoXML.org/2011/}geneRef" :
                generef_els.append(gene_etree)
        return generef_els

    def get_gene_via_generef(self, genesref_el):
        gene_els = []
        for generef_el in genesref_el:
            gene = fa.OrthoXMLQuery.getGeneFromId(generef_el.get("id"), self.Etree_root)
            if gene != None:
                gene_els.append(gene)
        return gene_els

    def write_OG(self, Etree_OG):
        if not self.cache_dir:
            raise RuntimeError('Check calling order.')

        # Create file var
        hog_nr = int(Etree_OG.get("id"))
        hog_id = "HOG{:06d}.orthoxml".format(hog_nr)
        fname = os.path.join(self.cache_dir, hog_id)
        print("Processing: ", fname)

        # Get element to store
        generef_els = self.get_generef_OG(Etree_OG)
        gene_els = self.get_gene_via_generef(generef_els)
        OG_el = Etree_OG

        # Get all information to store
        zoo = {} # <- {key:sp_etree || value: {key:db_el || values:[list_genes]}}
        for gene_el in gene_els: # <- for all gene el
            db_el = gene_el.getparent().getparent()
            sp_el = db_el.getparent()
            if sp_el in zoo.keys(): # <- if species already visited
                if db_el in zoo[sp_el].keys(): # <- if db already visited so add gene
                    zoo[sp_el][db_el].append(gene_el)
                else: # <- if db not visited so add db,genes
                    zoo[sp_el][db_el] = []
                    zoo[sp_el][db_el].append(gene_el)
            else: # <- if species not visited so add sp,db,gene
                zoo[sp_el] = {}
                zoo[sp_el][db_el] = []
                zoo[sp_el][db_el].append(gene_el)


        etree_2_dump = etree.Element("orthoXML")
        etree_2_dump.set("originVersion", 'Sep 2014')
        etree_2_dump.set("origin", 'OMA')
        etree_2_dump.set("version", '0.3')
        etree_2_dump.set("xmlns", 'http://orthoXML.org/2011/')


        for species_el in zoo.keys():
            species_xml = etree.Element("species")
            species_xml.set("name", species_el.get("name"))
            species_xml.set("NCBITaxId", species_el.get("NCBITaxId"))
            etree_2_dump.insert(0, species_xml)

            for db_el in zoo[species_el].keys():
                # Add <database> into <species>
                database_xml = etree.SubElement(species_xml, "database")
                database_xml.set("name", db_el.get("name"))
                database_xml.set("version", db_el.get("version"))

                # Add <genes> TAG into <database>
                genes_xml = etree.SubElement(database_xml, "genes")

                # Fill <genes> with <gene>
                for gene_el in zoo[species_el][db_el]:
                    gene_xml = etree.SubElement(genes_xml, "gene")
                    for attr, value in gene_el.attrib.items():
                        gene_xml.set(attr, value)

        groupsxml = etree.SubElement(etree_2_dump, "groups")
        str_OG_xml = etree.tostring(Etree_OG)
        str_OG_xml = str_OG_xml.decode().replace('xmlns="http://orthoXML.org/2011/"', "")
        OG_xml = etree.fromstring(str_OG_xml)
        groupsxml.append(OG_xml)

        indent(etree_2_dump)
        tree = etree.ElementTree(etree_2_dump)
        tree.write(fname, xml_declaration=True, encoding='utf-8', method="xml")
