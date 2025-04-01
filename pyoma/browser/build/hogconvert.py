import abc
import bisect
import collections
import copy
import logging
import os
import re
import string
from typing import Optional, Union

import numpy
import tables
import lxml.etree as etree
from ete3 import TreeNode

from .. import tablefmt
from ..convert import create_index_for_columns

HOGID_RE = re.compile(r"HOG:(?P<rel>[A-Z])?(?P<fam>\d+)(?:\.(?P<subhog>[a-z0-9.]*))?(?:_(?P<taxid>-?\d+))?")
logger = logging.getLogger(__name__)
Gene = collections.namedtuple("Gene", ["id", "e_nr", "original_id", "oma_id", "genome_taxonId"])
Species = collections.namedtuple("Species", ["species_code", "xml_taxonId", "entry_offset", "max_enr", "xml_node"])


class TaxonomyLookupHelper:
    def __init__(self, tree: TreeNode):
        self.taxtree = tree
        self._taxonId_to_node = {n.xml_taxonId: n for n in tree.traverse()}
        self._taxid_to_node = {n.taxid: n for n in tree.traverse()}

    def get_node_from_taxonId(self, taxon_id):
        return self._taxonId_to_node[int(taxon_id)]

    def get_node_from_taxid(self, taxid):
        return self._taxid_to_node[int(taxid)]

    def get_mrca(self, taxon_ids):
        taxon_ids = set(taxon_ids)
        leaves = [self._taxonId_to_node[taxonId] for taxonId in taxon_ids]
        if len(leaves) == 1:
            return leaves[0]
        mrca = self.taxtree.get_common_ancestor(*leaves)
        return mrca

    def levels_between(self, child_node, parent_node):
        n = child_node.up
        while n != parent_node:
            yield n
            n = n.up

    def as_xml(self):
        def _traverseR(node, xml):
            n = etree.SubElement(xml, "taxon", {"id": str(node.taxid), "name": str(node.name)})
            for child in node.children:
                _traverseR(child, n)

        root = etree.Element("taxonomy")
        _traverseR(self.taxtree, root)
        return root


class GeneLookupHelper:
    def __init__(self, parser: "AbstractOrthoXMLParser"):
        self._orig_generef = {}
        self._new_generef = {}
        self._taxonId_to_species = None
        self._entry_offsets = None
        self._species = []
        self.parser = parser

    def add_species(self, species, genes):
        self._orig_generef.update({int(x.id): x for x in genes})
        self._new_generef.update({int(x.e_nr): x for x in genes})
        self._species.append(species)

    def get_gene_by_original_id(self, original_id):
        return self._orig_generef.get(int(original_id))

    def get_gene_by_new_id(self, new_id):
        return self._new_generef.get(int(new_id))

    def _sort_and_freeze(self):
        self._species.sort(key=lambda x: x.entry_offset)
        self._taxonId_to_species = {x.xml_taxonId: x for x in self._species}
        self._entry_offsets = [x.entry_offset for x in self._species]

    def search_species_by_entry_nr(self, e_nr):
        if self._entry_offsets is None:
            self._sort_and_freeze()
        pos = bisect.bisect_left(self._entry_offsets, e_nr) - 1
        return self._species[pos]

    def iter_species_with_genes_nodes(self, new_generef_ids=None):
        genes = (
            (self._new_generef[enr] for enr in sorted(new_generef_ids))
            if new_generef_ids is not None
            else sorted(self._new_generef.values(), key=lambda x: x.e_nr)
        )

        sp, sp_node, genes_node = None, None, None
        for gene in genes:
            if sp is None or gene.e_nr > sp.max_enr:
                if sp is not None:
                    yield sp_node
                sp = self.search_species_by_entry_nr(gene.e_nr)
                sp_node = copy.deepcopy(sp.xml_node)
                sp_node.set("taxonId", str(self.parser.taxonomy.get_node_from_taxonId(sp.xml_taxonId).taxid))
                genes_node = sp_node.find(".//genes")

            assert sp.entry_offset < gene.e_nr <= sp.max_enr
            etree.SubElement(genes_node, "gene", {"id": str(gene.e_nr), "protId": gene.oma_id})
        yield sp_node


class Annotator:
    def __init__(self, parser: "AbstractOrthoXMLParser"):
        self.parser = parser
        self.parser.register_group_annotator(self)
        self._id_format = "HOG:{:s}{:07d}{:s}_{:s}"

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

    def update_hogids(self, hog: etree.Element, already_updated: bool = False, add_extra=False):
        if already_updated:
            taxon_accessor_func = self.parser.taxonomy.get_node_from_taxid
            gene_accessor_func = self.parser.gene_helper.get_gene_by_new_id
        else:
            taxon_accessor_func = self.parser.taxonomy.get_node_from_taxonId
            gene_accessor_func = self.parser.gene_helper.get_gene_by_original_id
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
                taxonId = node.get("taxonId")
                taxtree_node = taxon_accessor_func(taxonId)
                node.set("taxonId", str(taxtree_node.taxid))

                taxrange = node.find('./property[@name="TaxRange"]')
                if taxrange.get("value") != taxtree_node.name:
                    logger.warning("updating TaxRange property from %s to %s", taxrange.get("value"), taxtree_node.name)
                    taxrange.set("value", taxtree_node.name)

                orig_id = node.get("id")
                if orig_id is None or orig_id.isdigit() or (match := HOGID_RE.match(orig_id)) is None:
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
                    score_node.set("value", f"{len(covered_species) / taxtree_node.size:.3f}")
                else:
                    node.insert(
                        0,
                        etree.Element(
                            "score",
                            {"id": "CompletenessScore", "value": f"{len(covered_species) / taxtree_node.size:.3f}"},
                        ),
                    )
            elif node.tag == "paralogGroup":
                idx += 1
                next_og = f"{og}.{_getNextSubId(idx)}"
                for i, child in enumerate(node):
                    covered_species.update(_annotateGroupR(child, _encodeParalogClusterId(next_og, i), idx))
            elif node.tag == "geneRef":
                gene = gene_accessor_func(node.get("id"))
                node.set("id", str(gene.e_nr))
                if add_extra:
                    node.set("LOFT", og)
                covered_species.add(gene.genome_taxonId)
            return covered_species

        id_ = hog.get("id")
        if id_.isdigit():
            og = self.format_hogid(fam=id_)
        else:
            m = HOGID_RE.match(id_)
            if m is not None:
                og = self.format_hogid(m.group("fam"))
            else:
                raise ValueError(f"cannot parse roothog id of {hog}: {id_}")
        _annotateGroupR(hog, og, 0)
        return hog

    def annotate_hog(self, node):
        node = self.update_hogids(node, already_updated=False)
        return node

    def insert_missing_levels(self, node: etree._Element, tax_helper: TaxonomyLookupHelper):
        def insert_ogs_between(parent, child, tax_recent, tax_parent, nr_genes, nr_covered_species, include_self=True):
            pos = parent.index(child)
            if include_self:
                child = _insert_one_og(child, tax_recent, nr_covered_species, nr_genes)
            for lev in tax_helper.levels_between(tax_recent, tax_parent):
                child = _insert_one_og(child, lev, nr_covered_species=nr_covered_species, nr_genes=nr_genes)
            parent.insert(pos, child)

        def _insert_one_og(child, taxnode, nr_covered_species, nr_genes):
            el = etree.Element("orthologGroup", {"taxonId": str(taxnode.taxid)})
            el.append(
                etree.Element("score", {"id": "CompletenessScore", "value": f"{nr_covered_species / taxnode.size:.3f}"})
            )
            el.append(etree.Element("property", {"name": "TaxRange", "value": taxnode.name}))
            el.append(etree.Element("property", {"name": "taxid", "value": str(taxnode.taxid)}))
            el.append(etree.Element("property", {"name": "NrMemberGenes", "value": str(nr_genes)}))
            el.append(child)
            return el

        def recurse_traverse(node):
            if node.tag not in ("orthologGroup", "paralogGroup", "geneRef"):
                return
            if node.tag == "geneRef":
                genes = [self.parser.gene_helper.get_gene_by_new_id(node.get("id"))]
            else:
                genes = [self.parser.gene_helper.get_gene_by_new_id(n.get("id")) for n in node.findall(".//geneRef")]
            if node.tag != "orthologGroup":
                mrca = tax_helper.get_mrca([g.genome_taxonId for g in genes])
            else:
                mrca = tax_helper.get_node_from_taxid(node.get("taxonId"))
                for pos, c in enumerate(node.iterchildren()):
                    if c.tag not in ("score", "property"):
                        break
                node.insert(pos, etree.Element("property", {"name": "NrMemberGenes", "value": str(len(genes))}))
            try:
                parent = next(node.iterancestors("orthologGroup"))
                parent_tax_node = self.parser.taxonomy.get_node_from_taxid(parent.get("taxonId"))
            except StopIteration:
                parent_tax_node = None

            if parent_tax_node is not None:
                # insert missing levels as fake orthologGroup elements with the relevant data

                # Ortholog/Paralog Node - append missing tax range(s) as property tags under the current node
                if node.tag in ("orthologGroup", "paralogGroup"):
                    insert_ogs_between(
                        node.getparent(),
                        node,
                        mrca,
                        parent_tax_node,
                        len(genes),
                        len(set(g.genome_taxonId for g in genes)),
                        include_self=False,
                    )
                # GeneRef Node - insert ortholog node between self and parent; add all tax range(s) to new parent
                elif node.tag == "geneRef":
                    insert_ogs_between(node.getparent(), node, mrca, parent_tax_node, 1, 1, include_self=True)
                    return
            for child in node:
                recurse_traverse(child)

        recurse_traverse(node)
        return node

    def augment_hog(self, node):
        self.insert_missing_levels(node, self.parser.taxonomy)
        self.update_hogids(node, add_extra=True, already_updated=True)
        return node


class AbstractOrthoXMLParser(metaclass=abc.ABCMeta):
    def __init__(self, h5path):
        with tables.open_file(h5path, mode="r") as h5:
            self._read_from_h5(h5)
        self._taxname2taxid = {row["Name"].decode(): int(row["NCBITaxonId"]) for row in self._taxtab}
        self._orthoxml_attribs = {}
        self._scores = None
        self._group_observers = []
        self._group_annotators = []
        self.taxonomy = None
        self.gene_helper = GeneLookupHelper(self)

    def _read_from_h5(self, h5):
        self._taxtab = h5.get_node("/Taxonomy")[:]
        self._gstab = h5.get_node("/Genome")[:]
        self._rel_char = h5.get_node_attr("/", "oma_release_char")
        self.nr_proteins = len(h5.get_node("/Protein/Entries"))

    @property
    def rel_char(self):
        return self._rel_char

    def get_orthoxml_attribs(self):
        return self._orthoxml_attribs

    def register_group_handler(self, handler):
        self._group_observers.append(handler)

    def register_group_annotator(self, a):
        self._group_annotators.append(a)

    def signal_end_of_head(self):
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

        for annotator in self._group_annotators:
            node = annotator.augment_hog(node)

        for handler in self._group_observers:
            handler.process_augmented_hog(node)

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
        for n in root.traverse(strategy="preorder"):
            n.add_feature("size", len(n))
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

    def finish(self):
        for observer in self._group_observers:
            observer.finished()
        logger.debug("send finish signal to all observers")


class OrthoXMLGeneIdParserWithOmaProtId(AbstractOrthoXMLParser):
    def __init__(self, h5path: str):
        super().__init__(h5path)
        self._spcode2gs = {gs["UniProtSpeciesCode"].decode(): gs for gs in self._gstab}

    def process_species(self, node):
        assert node.tag == "{http://orthoXML.org/2011/}species"
        sp, sp_code, max_enr = None, None, 0
        xml_taxonId = int(node.get("taxonId"))
        genes = []
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
            genes.append(Gene(id_, e_nr, prot_id, prot_id, xml_taxonId))
        for genes_node in node.findall(".//{http://orthoXML.org/2011/}genes"):
            genes_node.clear()
        species = Species(sp_code, xml_taxonId, sp["EntryOff"], max_enr, strip_namespace(node))
        self.gene_helper.add_species(species, genes)
        logger.info(f"added {i + 1} genes for species {sp_code}")

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

    def _search_entry_nrs_of_ids_within_species(self, ids: numpy.array, sp: numpy.ndarray):
        en_min, en_max = sp["EntryOff"], sp["EntryOff"] + sp["TotEntries"]
        idx = numpy.searchsorted(self._prot_ids, ids, side="left", sorter=self._protkey)
        idy = numpy.searchsorted(self._prot_ids, ids, side="right", sorter=self._protkey)
        res = numpy.zeros(len(ids), dtype=int)
        for i, (low, high) in enumerate(zip(idx, idy)):
            sec = self._protkey[low:high]
            mask = (sec >= en_min) & (sec < en_max)
            match = numpy.where(mask)[0]
            if match.size != 1:
                logger.error(f"ID '{ids[i]}' not found or is not unique in the entry range [{en_min}:{en_max}]")
                logger.error(f" -> matches to {sec[match]}")
                raise ValueError(f"ID {ids[i]} not found or is not unique in the entry range [{en_min}:{en_max}]")
            res[i] = sec[match[0]]
        # entry numbers are 1-based
        return res + 1

    def process_species(self, node):
        assert node.tag == "{http://orthoXML.org/2011/}species"
        sp_name = node.get("name").encode("utf-8")
        try:
            sp = self._gstab[self._gstab["SciName"] == sp_name][0]
        except IndexError:
            try:
                sp = self._gstab[self._gstab["UniProtSpeciesCode"] == sp_name][0]
            except IndexError:
                logger.error(
                    f"no species found for {sp_name.decode()} in species table using SciName nor UniProtSpeciesCode"
                )
                raise RuntimeError(f"species {sp_name.decode()} not found.")
        sp_code = sp["UniProtSpeciesCode"].decode()
        xml_taxonId = int(node.get("taxonId"))

        genes = []
        sp_genes = {
            gene.get("id"): gene.get("protId").encode("utf-8")
            for gene in node.findall(".//{http://orthoXML.org/2011/}gene")
        }
        max_prot_id_len = min(max(len(z) for z in sp_genes.values()), max(len(z) for z in self._prot_ids))
        prot_ids = numpy.fromiter(sp_genes.values(), dtype=f"S{max_prot_id_len}")
        try:
            enrs = self._search_entry_nrs_of_ids_within_species(prot_ids, sp)
        except ValueError:
            logger.exception(
                f"Couldn't map {len(prot_ids)} protId from {sp_name.decode()} unambiguously to entry numbers"
            )
            raise

        genes = []
        for id_, prot_id, e_nr in zip(sp_genes.keys(), sp_genes.values(), enrs):
            oma_id = f"{sp_code}{e_nr-sp['EntryOff']:05d}"
            genes.append(Gene(id_, e_nr, prot_id, oma_id, xml_taxonId))

        for genes_node in node.findall(".//{http://orthoXML.org/2011/}genes"):
            genes_node.clear()
        species = Species(
            sp_code, xml_taxonId, sp["EntryOff"], sp["EntryOff"] + sp["TotEntries"], strip_namespace(node)
        )
        self.gene_helper.add_species(species, genes)
        logger.info(f"added {len(enrs)} genes for species {sp_code}")


def get_orthoxml_parser(h5path: Union[str, os.PathLike], is_oma_protId: bool = False) -> AbstractOrthoXMLParser:
    if is_oma_protId:
        return OrthoXMLGeneIdParserWithOmaProtId(h5path)
    else:
        return OrthoXMLGeneIdParserGeneralProtID(h5path)


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
        for sp in self._parser.gene_helper.iter_species_with_genes_nodes():
            self.writer.send(sp)
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


class FullAugmentedOrthoXMLObserver(FullOrthoXMLObserver):
    def process_hog(self, node):
        pass

    def process_augmented_hog(self, node: etree.Element):
        self.writer.send(node)


class PerFamilyHOGObserver(HogObserver):
    def __init__(
        self, parser: AbstractOrthoXMLParser, callback: Optional[callable], callback_augmented=Optional[callable]
    ):
        super().__init__(parser)
        self._callback = callback
        self._callback_augmented = callback_augmented

    def orthoxml_header_processed_hook(self):
        pass

    def _build_orthoxml(self, node):
        root_attribs = {"xmlns": "http://orthoXML.org/2011/", **self._parser.get_orthoxml_attribs()}
        root = etree.Element("orthoXML", attrib=root_attribs)
        genes = [int(n.get("id")) for n in node.findall(".//geneRef")]
        for sp in self._parser.gene_helper.iter_species_with_genes_nodes(genes):
            root.append(sp)
        root.append(self._parser.taxonomy.as_xml())
        root.append(self._parser.get_score_defs())
        groups = etree.SubElement(root, "groups")
        groups.append(node)
        orthoxml = etree.tostring(root, pretty_print=True, encoding="utf-8")
        return orthoxml

    def process_hog(self, node: etree.Element):
        if self._callback is not None:
            orthoxml = self._build_orthoxml(node)
            self._callback(node.get("id"), orthoxml)

    def process_augmented_hog(self, node: etree.Element):
        if self._callback_augmented is not None:
            orthoxml = self._build_orthoxml(node)
            self._callback_augmented(node.get("id"), orthoxml)


class HOGtoHDF5(HogObserver):
    def __init__(self, parser: AbstractOrthoXMLParser, h5path: str):
        super().__init__(parser)
        self.h5path = h5path
        self.nr_entries = parser.nr_proteins

    def __enter__(self):
        self.h5 = tables.open_file(self.h5path, mode="w", filters=tables.Filters(complib="blosc", complevel=6))
        orthoxml_group = self.h5.create_group("/", "OrthoXML")
        self.orthoxml_buffer = self.h5.create_earray(
            orthoxml_group, "Buffer", tables.StringAtom(1), (0,), "concatenated orthoxml files", expectedrows=1e9
        )
        self.orthoxml_buffer_augmented = self.h5.create_earray(
            orthoxml_group,
            "BufferAugmented",
            tables.StringAtom(1),
            (0,),
            "concatenated augmented orthoxml files",
            expectedrows=1e9,
        )
        self.orthoxml_index = self.h5.create_table(
            orthoxml_group,
            "Index",
            tablefmt.OrthoXmlHogTable,
            "Range index per HOG into OrthoXML Buffer",
            expectedrows=5e6,
        )
        self.leveltab = self.h5.create_table(
            "/",
            "HogLevel",
            tablefmt.HOGsTable,
            "nesting structure for each HOG",
            expectedrows=1e8,
        )
        self.index = {}
        self.hogid = numpy.zeros(shape=(self.nr_entries,), dtype=tablefmt.ProteinTable.columns["OmaHOG"].dtype)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.orthoxml_index.append(numpy.stack(list(self.index.values())))
        self.h5.create_carray("/", "OmaHOG", obj=self.hogid)
        self.h5.flush()
        if exc_type is None:
            # no exception happend. we build index and also per_level_tables
            create_index_for_columns(
                self.leveltab, "Fam", "ID", "Level", "CompletenessScore", "NrMemberGenes", "IsRoot"
            )
        self.h5.close()

    def process_augmented_hog(self, node: etree.Element):
        def get_hog_id(node):
            return node.get("id").split("_")[0]

        def get_hog_scores(og_node, tax_node):
            """extract the scores associated with an orthologGroup node

            only scores that are defined in HOGsTable are extract. The method
            returns a tuple with the scores in the order of the score fields."""
            all_score_ids = ("CompletenessScore", "ImpliedLosses")
            parse_fun = {"CompletenessScore": float, "ImpliedLosses": int}
            scores = collections.OrderedDict(
                [(score, tablefmt.HOGsTable.columns[score].dflt) for score in all_score_ids]
            )
            for score in og_node.iterfind("score"):
                score_id = score.get("id")
                scores[score_id] = parse_fun[score_id](score.get("value"))
            # might be overwritten in tax_node (more specific if available)
            for score_id in all_score_ids:
                val = tax_node.get(score_id)
                if val is not None:
                    scores[score_id] = parse_fun[score_id](val)
            return tuple(scores.values())

        def get_nr_member_genes(og_node):
            for child in og_node.iterfind("./property[@name='NrMemberGenes']"):
                return int(child.get("value"))
            logger.warning("couldn't find NrMemberGenes property. scanning xml file for geneRefs")
            return len(og_node.findall(".//geneRef"))

        match = HOGID_RE.match(node.get("id"))
        fam_nr = int(match.group("fam"))
        levs = []
        for taxnode in node.iterfind('.//property[@name="TaxRange"]'):
            ognode = taxnode.getparent()
            levs.append(
                (fam_nr, get_hog_id(ognode), taxnode.get("value"))
                + get_hog_scores(ognode, taxnode)
                + (
                    get_nr_member_genes(ognode),
                    bool(ognode.getparent().tag in ("paralogGroup", "groups")),
                    -1,  # default value for per taxlevel table index
                )
            )
        self.leveltab.append(levs)
        for gene in node.iterfind(".//geneRef"):
            self.hogid[int(gene.attrib["id"]) - 1] = gene.attrib["LOFT"].encode("utf-8")

    def _store_orthoxml_data(self, hogid, orthoxml, kind=None):
        if kind is None:
            buf, l_off, l_len = self.orthoxml_buffer, "HogBufferOffset", "HogBufferLength"
        elif kind == "augmented":
            buf, l_off, l_len = self.orthoxml_buffer_augmented, "HogAugmentedBufferOffset", "HogAugmentedBufferLength"
        match = HOGID_RE.match(hogid)
        fam = int(match.group("fam"))
        off = len(buf)
        buf.append(numpy.ndarray((len(orthoxml),), buffer=orthoxml, dtype=tables.StringAtom(1)))
        if fam not in self.index:
            self.index[fam] = numpy.zeros(1, dtype=self.orthoxml_index.dtype)
        e = self.index[fam]
        e["Fam"] = fam
        e[l_off] = off
        e[l_len] = len(orthoxml)

    def store_orthoxml(self, hogid, orthoxml):
        self._store_orthoxml_data(hogid, orthoxml)

    def store_orthoxml_augmented(self, hogid, orthoxml):
        self._store_orthoxml_data(hogid, orthoxml, kind="augmented")


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
            elif elem.tag == "{http://orthoXML.org/2011/}scores":
                handler.process_scores(elem)
            elif elem.tag == "{http://orthoXML.org/2011/}orthoXML":
                logger.info(f"Finished parsing orthoxml file. processed {nr_roothogs_parsed} roothogs.")
    handler.finish()
