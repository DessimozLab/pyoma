import familyanalyzer as fa
import lxml.etree as etree
import os


class OrthoXMLSplitter(object):
    def __init__(self, xml_file):
        self.xml_file = xml_file
        self.Etree_XML = etree.parse(self.xml_file)
        self.Etree_root = self.Etree_XML.getroot()
        self.Etree_OGs = fa.OrthoXMLQuery.getToplevelOrthologGroups(self.Etree_root)
        self.gene_lookup = {gene.get('id'): gene for gene in fa.OrthoXMLQuery.getInputGenes(self.Etree_root)}

    def split_each_hog_into_individual(self, storage_folder, list_hog_nr=None):
        os.system("mkdir " + storage_folder)
        os.system("rm " + storage_folder + '/*')
        for og in self.Etree_OGs:
            hog_nr = int(og.get("id"))
            hog_id = "HOG{:06d}.orthoxml".format(hog_nr)
            fname = os.path.join(storage_folder, hog_id)
            print("Processing: ", fname)

            if list_hog_nr is not None:
                if int(og.get("id")) in list_hog_nr:
                    self.create_new_orthoxml(fname, [og])
            else:
                self.create_new_orthoxml(fname, [og])

    def extract_hogs_into_new_orthoxml(self, new_fn, list_hog_nr):
        ogs = []
        for og_etree in self.Etree_OGs:
            if str(og_etree.get("id")) in list_hog_nr:
                ogs.append(og_etree)
        self.create_new_orthoxml(new_fn, ogs)

    def get_generef_OG(self, og_node):
        return fa.OrthoXMLQuery.getGeneRefNodes(og_node, recursively=True)

    def get_gene_via_generef(self, genesref_ids):
        genesref_ids = set(genesref_ids)
        return [self.gene_lookup[gene_id] for gene_id in genesref_ids]

    def create_new_orthoxml(self, fn, OGs):
        # Get element to store
        for og_node in OGs:
            gene_ids = [gene_ref_elem.get("id") for gene_ref_elem in self.get_generef_OG(og_node)]
        gene_els = self.get_gene_via_generef(gene_ids)

        # Get all information to store
        zoo = {}  # <- {key:sp_etree || value: {key:db_el || values:[list_genes]}}
        for gene_el in gene_els:  # <- for all gene el
            db_el = gene_el.getparent().getparent()
            sp_el = db_el.getparent()
            if sp_el in zoo.keys():  # <- if species already visited
                if db_el in zoo[sp_el].keys():  # <- if db already visited so add gene
                    zoo[sp_el][db_el].append(gene_el)
                else:  # <- if db not visited so add db,genes
                    zoo[sp_el][db_el] = []
                    zoo[sp_el][db_el].append(gene_el)
            else:  # <- if species not visited so add sp,db,gene
                zoo[sp_el] = {}
                zoo[sp_el][db_el] = []
                zoo[sp_el][db_el].append(gene_el)

        etree_2_dump = etree.Element("orthoXML")
        for attr, value in self.Etree_root.items():
            etree_2_dump.set(attr, value)

        for species_el in zoo.keys():
            species_xml = etree.Element("species")
            for attr, value in species_el.items():
                species_xml.set(attr, value)
            etree_2_dump.insert(0, species_xml)

            for db_el in zoo[species_el].keys():
                # Add <database> into <species>
                database_xml = etree.SubElement(species_xml, "database")
                for attr, value in db_el.items():
                    database_xml.set(attr, value)

                # Add <genes> TAG into <database>
                genes_xml = etree.SubElement(database_xml, "genes")

                # Fill <genes> with <gene>
                for gene_el in zoo[species_el][db_el]:
                    gene_xml = etree.SubElement(genes_xml, "gene")
                    for attr, value in gene_el.attrib.items():
                        gene_xml.set(attr, value)

        groupsxml = etree.SubElement(etree_2_dump, "groups")
        for og_et in OGs:
            str_OG_xml = etree.tostring(og_et)
            str_OG_xml = str_OG_xml.decode().replace('xmlns="http://orthoXML.org/2011/"', "")
            OG_xml = etree.fromstring(str_OG_xml)
            groupsxml.append(OG_xml)

        tree = etree.ElementTree(etree_2_dump)
        tree.write(fn, xml_declaration=True, encoding='utf-8', method="xml", pretty_print=True)

