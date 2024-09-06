import abc
import collections
import copy
import logging
import re
import string
from typing import Optional

import tables
import lxml.etree as etree
from ete3 import TreeNode

logger = logging.getLogger(__name__)
Gene = collections.namedtuple("Gene", ["id", "e_nr", "original_id", "oma_id", "genome_taxonId"])
Species = collections.namedtuple("Species", ["species_code", "xml_taxonId", "entry_offset", "max_enr", "xml_node"])


class TaxonomyLookupHelper:
    def __init__(self, tree: TreeNode):
        self.taxtree = tree
        self._taxonId_to_node = {n.xml_taxonId: n for n in tree.traverse()}

    def get_node_from_taxonId(self, taxon_id):
        return self._taxonId_to_node[int(taxon_id)]

    def as_xml(self):
        def _traverseR(node, xml):
            n = etree.SubElement(xml, "taxon", {"id": str(node.taxid), "name": str(node.name)})
            for child in node.children:
                _traverseR(child, n)

        root = etree.Element("taxonomy")
        _traverseR(self.taxtree, root)
        return root


class Annotator:
    def __init__(self, parser):
        self.parser = parser
        self.parser.register_group_annotator(self)
        self._id_format = "HOG:{:s}{:07d}{:s}_{:s}"
        self._hogid_re = re.compile(
            r"HOG:(?P<rel>[A-Z])?(?P<fam>\d+)(?:\.(?P<subhog>[a-z0-9.]*))?(?:_(?P<taxid>-?\d+))?"
        )

    def format_hogid(self, fam, subhog=None, taxid=None):
        num = f"{int(fam):07d}"
        id_ = f"HOG:{self.parser.rel_char}{num}"
        if subhog is not None:
            id_ += f".{subhog}"
        if taxid is not None:
            id_ += f"_{taxid}"
        return id_

    def provide_scores(self):
        return [
            ("CompletenessScore", "Fraction of expected species genes in the HOG"),
        ]

    def update_hogids(self, hog: etree.Element):
        dupCnt = collections.defaultdict(int)

        def _getNextSubId(idx):
            """helper method to return the next number at a given depth of
            duplication (idx)"""
            dupCnt[idx - 1] += 1
            return dupCnt[idx - 1]

        def _encodeParalogClusterId(prefix, nr):
            letters = []
            while nr // 26 > 0:
                letters.append(string.ascii_lowercase[nr % 26])
                nr = nr // 26 - 1
            letters.append(string.ascii_lowercase[nr % 26])
            return prefix + "".join(letters[::-1])

        def _annotateGroupR(node: etree.Element, og: str, idx: int = 0) -> set:
            """create the og attributes at the orthologGroup elements
            according to the naming schema of LOFT. ParalogGroup elements
            do not get own attributes (not possible in the xml schema),
            but propagate their sub-names for the subsequent orthologGroup
            elements."""
            covered_species = set([])
            if node.tag == "orthologGroup":
                orig_xml_taxonId = node.get("taxonId")
                taxtree_node = self.parser.taxonomy.get_node_from_taxonId(orig_xml_taxonId)
                node.set("taxonId", str(taxtree_node.taxid))

                taxrange = node.find('./property[@name="TaxRange"]')
                if taxrange.get("value") != taxtree_node.name:
                    logger.info("updating TaxRange property from %s to %s", taxrange.get("value"), taxtree_node.name)
                    taxrange.set("value", taxtree_node.name)

                orig_id = node.get("id")
                if orig_id is None or orig_id.isdigit() or (match := self._hogid_re.match(orig_id) is None):
                    node.set("id", f"{og}_{taxtree_node.taxid}")
                else:
                    id_ = self.format_hogid(match.group("fam"), match.group("subhog"), taxtree_node.taxid)
                    if og != id_.split("_", maxsplit=1)[0]:
                        logger.warning("expected subfamily id does not match existing ID.")
                    node.set("id", id_)

                for child in node:
                    covered_species.update(_annotateGroupR(child, og, idx))
                score_node = node.find('./score[@id="CompletenessScore"]')
                if score_node is not None:
                    score_node.set("value", f"{len(covered_species) / len(taxtree_node):.3f}")
                else:
                    node.insert(
                        0,
                        etree.Element(
                            "score",
                            {"id": "CompletenessScore", "value": f"{len(covered_species) / len(taxtree_node):.3f}"},
                        ),
                    )
            elif node.tag == "paralogGroup":
                idx += 1
                next_og = f"{og}.{_getNextSubId(idx)}"
                for i, child in enumerate(list(node)):
                    covered_species.update(_annotateGroupR(child, _encodeParalogClusterId(next_og, i), idx))
            elif node.tag == "geneRef":
                gene = self.parser.get_gene(node.get("id"))
                node.set("id", str(gene.e_nr))
                covered_species.add(gene.genome_taxonId)
            return covered_species

        id_ = hog.get("id")
        if id_.isdigit():
            og = self.format_hogid(fam=id_)
        else:
            m = self._hogid_re.match(id_)
            if m is not None:
                og = self.format_hogid(m.group("fam"))
            else:
                raise ValueError(f"cannot parse roothog id of {hog}: {id_}")
        _annotateGroupR(hog, og, 0)
        return hog

    def annotate_hog(self, node):
        logger.debug("before annotate_hog: " + etree.tostring(node).decode("utf-8"))
        node = self.update_hogids(node)
        logger.debug("after annotate_hog: " + etree.tostring(node).decode("utf-8"))
        return node


class AbstractOrthoXMLParser(metaclass=abc.ABCMeta):
    def __init__(self, h5path):
        with tables.open_file(h5path, mode="r") as h5:
            self._read_from_h5(h5)
        self._taxname2taxid = {row["Name"].decode(): int(row["NCBITaxonId"]) for row in self._taxtab}
        self._taxonomyId2taxid = {}
        self._geneid2gene = {}
        self._orthoxml_attribs = {}
        self._taxonId2species = {}
        self._scores = None
        self._group_observers = []
        self._group_annotators = []
        self.taxonomy = None

    def _read_from_h5(self, h5):
        self._taxtab = h5.get_node("/Taxonomy")[:]
        self._gstab = h5.get_node("/Genome")[:]
        self._rel_char = h5.get_node_attr("/", "oma_release_char")

    @property
    def rel_char(self):
        return self._rel_char

    def get_orthoxml_attribs(self):
        return self._orthoxml_attribs

    # def get_species_node(self, xml_taxonId=None, taxid=None):
    #     sp = copy.deepcopy(self._taxonId2species[str(xml_taxonId)])
    #     sp.xml_node.set('taxonId', str(self.taxonomy.get_node_from_taxonId(sp.xml_taxonId).taxid))

    def iter_species_nodes(self):
        for sp in self._taxonId2species.values():
            sp = copy.deepcopy(sp)
            sp.xml_node.set("taxonId", str(self.taxonomy.get_node_from_taxonId(sp.xml_taxonId).taxid))
            yield sp

    def iter_genes(self):
        yield from self._geneid2gene.values()

    def get_gene(self, generef_id):
        return self._geneid2gene[int(generef_id)]

    def register_group_handler(self, handler):
        self._group_observers.append(handler)

    def register_group_annotator(self, a):
        self._group_annotators.append(a)

    def signal_end_of_head(self):
        self._taxonId2species = {
            k: v for (k, v) in sorted(self._taxonId2species.items(), key=lambda x: x[1].entry_offset)
        }
        self._geneid2gene = {k: v for (k, v) in sorted(self._geneid2gene.items(), key=lambda x: x[1].e_nr)}
        for handler in self._group_observers:
            handler.orthoxml_header_processed_hook()

    def process_orthoxml_head(self, node):
        self._orthoxml_attribs = node.attrib

    @abc.abstractmethod
    def process_species(self, node):
        pass

    def process_group(self, node):
        assert node.tag == "{http://orthoXML.org/2011/}orthologGroup"
        node = strip_namespace(node)
        for annotator in self._group_annotators:
            node = annotator.annotate_hog(node)

        for handler in self._group_observers:
            handler.process_hog(node)

    def _get_features_for_taxonomy(self, xml_node: etree.Element):
        xml_taxonId = int(xml_node.get("id"))
        name = xml_node.get("name")
        taxid = self._taxname2taxid[name]
        return {"taxid": taxid, "name": name, "xml_taxonId": xml_taxonId}

    def _build_tree_node(self, xml_node: etree.Element, parent: Optional[TreeNode] = None) -> TreeNode:
        tn = TreeNode()
        tn.add_features(**self._get_features_for_taxonomy(xml_node))
        if parent is not None:
            parent.add_child(tn)
        return tn

    def _post_process_taxonomy_tree(self, tree_node: TreeNode):
        return tree_node

    def process_taxonomy(self, node):
        assert node.tag == "{http://orthoXML.org/2011/}taxonomy"

        def traverse(xml_node: etree.Element, tree: TreeNode):
            for child in xml_node.iterchildren(tag="{http://orthoXML.org/2011/}taxon"):
                child_node = self._build_tree_node(child, tree)
                traverse(child, child_node)

        root_tax = node.findall("./{http://orthoXML.org/2011/}taxon")
        assert len(root_tax) == 1
        root_tax = root_tax[0]
        root = self._build_tree_node(root_tax)
        traverse(root_tax, root)
        root = self._post_process_taxonomy_tree(root)
        self.taxonomy = TaxonomyLookupHelper(root)

    def process_scores(self, node):
        assert node.tag == "{http://orthoXML.org/2011/}scores"
        self._scores = strip_namespace(node)

    def get_score_defs(self):
        score_ids = set([])
        if self._scores is not None:
            score_ids = set(n.get("id") for n in self._scores.iter("scoreDef"))
        else:
            self._scores = etree.Element("scores")
        for annotator in self._group_annotators:
            for score_id, score_def in annotator.provide_scores():
                if score_id not in score_ids:
                    etree.SubElement(self._scores, "scoreDef", id=score_id, desc=score_def)
        return self._scores if len(self._scores) > 0 else None


class OrthoXMLGeneIdParserWithOmaProtId(AbstractOrthoXMLParser):
    def __init__(self, h5path: str):
        super().__init__(h5path)
        self._spcode2gs = {gs["UniProtSpeciesCode"].decode(): gs for gs in self._gstab}

    def process_species(self, node):
        assert node.tag == "{http://orthoXML.org/2011/}species"
        sp, sp_code, max_enr = None, None, 0
        xml_taxonId = int(node.get("taxonId"))
        for i, gene in enumerate(node.findall(".//{http://orthoXML.org/2011/}gene")):
            id_ = int(gene.get("id"))
            prot_id = gene.get("protId")
            if sp is None:
                sp_code, nr = prot_id[:5], int(prot_id[5:])
                sp = self._spcode2gs[sp_code]
            else:
                nr = int(gene.get("protId")[5:])
            e_nr = int(sp["EntryOff"]) + nr
            max_enr = max(max_enr, e_nr)
            self._geneid2gene[id_] = Gene(id_, e_nr, prot_id, prot_id, xml_taxonId)
        for genes_node in node.findall(".//{http://orthoXML.org/2011/}genes"):
            genes_node.clear()
        logger.info(f"added {i+1} genes for species {sp_code}")
        self._taxonId2species[xml_taxonId] = Species(sp_code, xml_taxonId, sp["EntryOff"], max_enr, node)

    def _get_features_for_taxonomy(self, xml_node: etree.Element):
        xml_taxonId = int(xml_node.get("id"))
        name = xml_node.get("name")
        try:
            taxid = self._taxname2taxid[name]
            code = None
        except KeyError:
            gs = self._spcode2gs[name]
            name = gs["SciName"].decode()
            taxid = int(gs["NCBITaxonId"])
            code = gs["UniProtSpeciesCode"].decode()
        res = {"taxid": taxid, "name": name, "xml_taxonId": xml_taxonId}
        if code is not None:
            res["code"] = code
        return res

    def _post_process_taxonomy_tree(self, tree_node):
        nxt_fix_taxid = -10
        for leaf in tree_node.iter_leaves():
            if (leaf.taxid == leaf.up.taxid) or (leaf.name == leaf.up.name):
                leaf.name = leaf.code
                leaf.taxid = nxt_fix_taxid
                nxt_fix_taxid -= 1
        taxid_counter = collections.Counter(l.taxid for l in tree_node.traverse())
        if max(taxid_counter.values()) > 1:
            raise ValueError("taxids are not unique: ", taxid_counter.most_common(3))
        return tree_node


class OrthoXMLGeneIdParserGeneralProtID(OrthoXMLGeneIdParserWithOmaProtId):
    def __init__(self, h5path):
        super().__init__(h5path)
        self._protkey = self._prot_ids.argsort()

    def _read_from_h5(self, h5):
        super()._read_from_h5(h5)
        self._prot_ids = h5.get_node("/Protein/Entries").read(field="CanonicalId")


class HogObserver:
    def __init__(self, parser: AbstractOrthoXMLParser):
        self._parser = parser
        parser.register_group_handler(self)

    def orthoxml_header_processed_hook(self):
        pass

    def process_hog(self, node: etree.Element):
        pass

    def process_augmented_hog(self, node: etree.Element):
        pass

    def finished(self):
        pass


def strip_namespace(element):
    """Remove namespace from the given element and its descendants."""
    new_element = etree.Element(etree.QName(element).localname, nsmap=None)
    # Copy attributes
    new_element.attrib.update(element.attrib)
    # Recursively copy children
    for child in element:
        new_child = strip_namespace(child)
        new_element.append(new_child)

    # If element has text, add it
    if element.text:
        new_element.text = element.text
    return new_element


def xml_writer(stream, root_elem):
    with open(stream, "wb") as xf:
        xf.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
        # Write the root element's opening tag
        xf.write(b"<%s" % root_elem["tag"].encode("utf-8"))
        if root_elem["xmlns"] is not None:
            xf.write(b' xmlns="%s"' % root_elem["xmlns"].encode("utf-8"))
        for key, value in root_elem["attribs"].items():
            xf.write(b' %s="%s"' % (key.encode("utf-8"), value.encode("utf-8")))
        xf.write(b">\n")

        while True:
            try:
                data = yield
            except GeneratorExit:
                break
            if isinstance(data, etree._Element):
                # Write an XML element
                xf.write(etree.tostring(strip_namespace(data), pretty_print=True, encoding="utf-8"))
            elif isinstance(data, (str, bytes)):
                # Write raw string/byte data
                if isinstance(data, str):
                    data = data.encode("utf-8")
                xf.write(data)
            else:
                raise ValueError("Unsupported data type sent to writer")
            xf.flush()
        # Close the root element
        xf.write(b"</%s>\n" % root_elem["tag"].encode("utf-8"))


class FullOrthoXMLObserver(HogObserver):
    def __init__(self, parser: AbstractOrthoXMLParser, outpath: str):
        super().__init__(parser)
        self._outpath = outpath
        self.writer = None

    def _write_xml_header(self):
        root_elem = {
            "tag": "orthoXML",
            "xmlns": "http://orthoXML.org/2011/",
            "attribs": self._parser.get_orthoxml_attribs(),
        }
        self.writer = xml_writer(self._outpath, root_elem)
        next(self.writer)

    def orthoxml_header_processed_hook(self):
        self._write_xml_header()
        sp_it = self._parser.iter_species_nodes()
        sp, genes = None, None
        for gene in self._parser.iter_genes():
            if sp is None or sp.max_enr < gene.e_nr:
                if sp is not None:
                    self.writer.send(sp_node)
                    sp_node = None
                sp = next(sp_it)
                sp_node = sp.xml_node
                genes = sp_node.find(".//{http://orthoXML.org/2011/}genes")
            etree.SubElement(genes, "gene", {"id": str(gene.e_nr), "protId": gene.oma_id})
        self.writer.send(sp_node)
        sp_node = None
        self.writer.send(self._parser.taxonomy.as_xml())
        scores = self._parser.get_score_defs()
        if scores is not None:
            self.writer.send(scores)
        self.writer.send("<groups>")

    def process_hog(self, node: etree.Element):
        self.writer.send(node)

    def finished(self):
        self.writer.send("</groups>")
        self.writer.close()
        super().finished()


def parse_orthoxml(xml, handler):
    og_depth = 0
    nr_roothogs_parsed = 0
    for event, elem in etree.iterparse(xml, events=("start", "end")):
        if event == "start":
            if elem.tag == "{http://orthoXML.org/2011/}orthologGroup":
                og_depth += 1
            elif elem.tag == "{http://orthoXML.org/2011/}orthoXML":
                handler.process_orthoxml_head(elem)
            elif elem.tag == "{http://orthoXML.org/2011/}groups":
                handler.signal_end_of_head()
        elif event == "end":
            if elem.tag == "{http://orthoXML.org/2011/}orthologGroup":
                og_depth -= 1
                if og_depth == 0:
                    handler.process_group(elem)
                    elem.clear()
                    nr_roothogs_parsed += 1
                    if nr_roothogs_parsed % 10000 == 0:
                        logger.info("parsed %d rootlevel groups", nr_roothogs_parsed)
            elif elem.tag == "{http://orthoXML.org/2011/}species":
                handler.process_species(elem)
            elif elem.tag == "{http://orthoXML.org/2011/}taxonomy":
                handler.process_taxonomy(elem)
                elem.clear()
            elif elem.tag == "{http://orthoXML.org/2011/}scoreDef":
                handler.process_score_def(elem)
