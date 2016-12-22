import familyanalyzer as fa
import lxml.etree as etree
import os


class OrthoXMLSplitter(object):
    def __init__(self, xml_file):
        self.xml_file = xml_file
        self.Etree_XML = etree.parse(self.xml_file)
        self.Etree_root = self.Etree_XML.getroot()
        self.Etree_OGs = fa.OrthoXMLQuery.getToplevelOrthologGroups(self.Etree_root)
        self.Etree_header_genes = fa.OrthoXMLQuery.getInputGenes(self.Etree_root)

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
                pass
                self.create_new_orthoxml(fname, [og])

    def extract_hogs_into_new_orthoxml(self, new_fn, list_hog_nr):
        ogs = []
        for og_etree in self.Etree_OGs:
            if str(og_etree.get("id")) in list_hog_nr:
                ogs.append(og_etree)
        self.create_new_orthoxml(new_fn, ogs)

    def get_generef_OG(self, Etree_OG):
        generef_els = []
        for gene_etree in Etree_OG.getiterator():
            if gene_etree.tag == "{http://orthoXML.org/2011/}geneRef":
                generef_els.append(gene_etree)
        return generef_els

    def get_gene_via_generef(self, genesref_ids):
        genesref_ids = set(genesref_ids)
        gene_els = []
        cpt = 0
        for gene_header in self.Etree_header_genes:
            if gene_header.get("id") in genesref_ids:
                gene_els.append(gene_header)
        return gene_els

    def indent(self, elem, level=0):
        """
        re structure the xml tree in human readable format (pre processing before writing the tree in a file)
        :param elem:
        :param level:
        :return:
        """
        i = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                self.indent(elem, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i

    def create_new_orthoxml(self, fn, OGs):

        # Get element to store
        generef_ids = []
        for etree_og in OGs:
            og_gnref = self.get_generef_OG(etree_og)
            for gnref in og_gnref:
                generef_ids.append(gnref.get("id"))
        gene_els = self.get_gene_via_generef(generef_ids)

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
        for og_et in OGs:
            str_OG_xml = etree.tostring(og_et)
            str_OG_xml = str_OG_xml.decode().replace('xmlns="http://orthoXML.org/2011/"', "")
            OG_xml = etree.fromstring(str_OG_xml)
            groupsxml.append(OG_xml)

        self.indent(etree_2_dump)
        tree = etree.ElementTree(etree_2_dump)
        tree.write(fn, xml_declaration=True, encoding='utf-8', method="xml")

    def write_OG_(self, Etree_OG):

        # Create file var
        hog_nr = int(Etree_OG.get("id"))
        hog_id = "HOG{:06d}.html".format(hog_nr)
        fname = os.path.join(self.dir_storage, hog_id)
        print("Processing: ", fname)

        # Get element to store
        generef_els = self.get_generef_OG(Etree_OG)
        gene_els = self.get_gene_via_generef(generef_els)
        OG_el = Etree_OG

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

        self.indent(etree_2_dump)
        tree = etree.ElementTree(etree_2_dump)
        tree.write(fname, xml_declaration=True, encoding='utf-8', method="xml")
