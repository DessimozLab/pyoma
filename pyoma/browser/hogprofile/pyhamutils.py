import pyham
import xml.etree.ElementTree as ET
import logging

logger = logging.getLogger(__name__)


def get_ham_treemap_from_row(row, tree):
    """get the mapping of a pyham tree profile for a fiven root hog"""
    fam, orthoxml = row
    orthoxml = switch_name_ncbi_id(orthoxml)
    try:
        ham_obj = pyham.Ham(
            tree,
            orthoxml,
            type_hog_file="orthoxml",
            use_internal_name=True,
            orthoXML_as_string=True,
        )
        tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0])
        return tp.treemap
    except (TypeError, AttributeError) as err:
        logger.exception("Error in pyham treeprofile occured")
        return None


def switch_name_ncbi_id(orthoxml):
    """swap ncbi taxid for species name to avoid ambiguity"""
    root = ET.fromstring(orthoxml)
    for node in root:
        if "species" in node.tag:
            node.attrib["name"] = node.attrib["NCBITaxId"]
    orthoxml = ET.tostring(root, encoding="unicode", method="xml")
    return orthoxml
