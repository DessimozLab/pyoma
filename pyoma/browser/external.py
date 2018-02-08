import collections
import os

import math
from tqdm import tqdm
from lxml import etree
from .db import Database


"""Module to generate external data and crosslinks for NCBI.

"""


class NCBILinkOutXML(object):
    root_node = "not_set"
    provider_id = "9822"

    def __init__(self):
        root = etree.Element(self.root_node)
        for key, value in self.root_children().items():
            root.append(self.text_elemement(key, value))
        self.tree = etree.ElementTree(root)
        self._add_doctype()

    def root_children(self):
        return {}

    def _add_doctype(self):
        self.tree.docinfo.public_id = '-//NLM//DTD LinkOut 1.0//EN'
        self.tree.docinfo.system_url = 'https://www.ncbi.nlm.nih.gov/projects/linkout/doc/LinkOut.dtd'

    def text_elemement(self, tag, text):
        el = etree.Element(tag)
        el.text = text
        return el

    def write(self, fh):
        fh.write(etree.tostring(self.tree, pretty_print=True, xml_declaration=True))


class Provider(NCBILinkOutXML):
    root_node = "Provider"

    def root_children(self):
        elements = collections.OrderedDict(
            [("ProviderId", self.provider_id),
             ("Name", "OMA Browser: Orthologous MAtrix"),
             ("NameAbbr", "OMA"),
             ("SubjectType", "taxonomy/phylogenetic"),
             ("Url", "http://omabrowser.org"),
             ("Brief", "OMA is a method and database for the inference of orthologs among complete genomes. "
                       "We provide browsable orthology predictions, APIs, flat file downloads among thousands "
                       "of genomes.")])
        return elements


class Resource(NCBILinkOutXML):
    root_node = "LinkSet"
    lnk_nr = 1

    def _add_objs(self, accs):
        objsel = etree.Element("ObjectSelector")
        objsel.append(self.text_elemement("Database", self.database()))
        objlst = etree.Element("ObjectList")
        objsel.append(objlst)
        for acc in accs:
            objlst.append(self.text_elemement("ObjId", acc))
        return objsel

    def _add_url_section(self, acc):
        el = etree.Element('ObjectUrl')
        el.append(self.text_elemement('Base', self.base_url()))
        nxt = rule = etree.Element("Rule")
        for k, rule_part in enumerate(self.rule_url(acc)):
            if isinstance(rule_part, str):
                if k == 0:
                    nxt.text = rule_part
                else:
                    nxt.tail = rule_part
            elif rule_part.tag == etree.Entity:
                nxt.append(rule_part)
                nxt = rule_part
        el.append(rule)
        el.append(self.text_elemement('SubjectType', "taxonomy/phylogenetic"))
        return el

    def add_link(self, accs):
        lnk = etree.Element("Link")
        lnk.append(self.text_elemement('LinkId', str(self.lnk_nr)))
        lnk.append(self.text_elemement('ProviderId', self.provider_id))
        lnk.append(self._add_objs(accs))
        lnk.append(self._add_url_section(accs))
        self.tree.getroot().append(lnk)
        self.lnk_nr += 1

    def database(self):
        return "not set"

    def base_url(self):
        return "https://omabrowser.org/oma/hogs/"

    def rule_url(self, acc):
        return "",


class GenesResource(Resource):
    DISKSIZE_HEADER = 200
    DISKSIZE_PER_LINK = 435
    base_name = 'resource_genes'

    def base_url(self):
        return "https://omabrowser.org/cgi-bin/gateway.pl/"

    def rule_url(self, acc):
        return "?f=DisplayEntry&p1=" + next(iter(acc.values())),

    def database(self):
        return "Gene"


class ProteinResource(Resource):
    DISKSIZE_HEADER = 500
    DISKSIZE_PER_LINK = 40
    base_name = 'resource_protein'

    def base_url(self):
        return "https://omabrowser.org/oma/hogs/"

    def rule_url(self, acc):
        return etree.Entity("lo.pacc"), "/vis/"

    def database(self):
        return "Protein"


class TaxonomyResource(Resource):
    DISKSIZE_HEADER = 200
    DISKSIZE_PER_LINK = 435
    base_name = 'resource_taxonomy'

    def database(self):
        return "taxonomy"

    def base_url(self):
        return "https://omabrowser.org/cgi-bin/gateway.pl/"

    def rule_url(self, acc):
        return "?f=DisplayOS&p1=" + next(iter(acc.values())),


class LinkoutBuffer(object):
    def __init__(self, resource, outdir, bulk_add=True, max_file_size=20*2**20):
        self.max_records = math.floor((max_file_size - resource.DISKSIZE_HEADER) /
                                      resource.DISKSIZE_PER_LINK)
        self.cur_nr = 0
        self.buf = []
        self.bulk_add = bulk_add
        self.resource_type = resource
        self.outdir = outdir

    def add(self, obj):
        self.buf.append(obj)
        if len(self.buf) >= self.max_records:
            self.flush()

    def flush(self):
        res = self.resource_type()
        if self.bulk_add:
            res.add_link(self.buf)
        else:
            for obj in self.buf:
                res.add_link(obj)
        fn = os.path.join(self.outdir,
                          '{}_{:03d}.xml'.format(res.base_name, self.cur_nr))
        with open(fn, 'wb') as fh:
            res.write(fh)
        self.cur_nr += 1
        self.buf = []


def run(outdir='/tmp', infile='../pyomabrowser/OmaServer.h5'):
    prov = Provider()
    with open(os.path.join(outdir, 'provider.xml'), 'wb') as fh:
        prov.write(fh)

    db = Database(infile)
    xrefs = db.get_hdf5_handle().get_node('/XRef')
    xref_source_enum = xrefs.get_enum('XRefSource')

    protein_buffer = LinkoutBuffer(ProteinResource, outdir, bulk_add=True)
    genes_buffer = LinkoutBuffer(GenesResource, outdir, bulk_add=False)
    for xref in tqdm(xrefs):
        if xref['XRefSource'] == xref_source_enum['RefSeq']:
            protein_buffer.add(xref['XRefId'].decode())
        elif xref['XRefSource'] == xref_source_enum['EntrezGene']:
            genes_buffer.add({xref['XRefId'].decode(): db.id_mapper['OMA'].map_entry_nr(xref['EntryNr'])})
    protein_buffer.flush()
    genes_buffer.flush()

    with open(os.path.join(outdir, 'resource_taxonomy.xml'), 'wb') as fh:
        taxs = TaxonomyResource()
        for row in db.id_mapper['OMA'].genome_table:
            taxs.add_link({str(row['NCBITaxonId']): row['UniProtSpeciesCode'].decode()})
        taxs.write(fh)
