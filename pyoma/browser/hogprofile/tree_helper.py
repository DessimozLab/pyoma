import copy
import collections
import itertools
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import ete4


def get_newick_tree_from_tax_db(tax):
    def traverse(node):
        yield node
        if "children" in node:
            for child in node["children"]:
                yield from traverse(child)

    def rec(node):
        if "children" not in node:
            return str(node["id"])
        children = []
        for child in node["children"]:
            children.append(rec(child))
        return "({}){}".format(",".join(children), node["id"])

    def get_duplicates(node):
        c = collections.Counter(x["id"] for x in traverse(node))
        return list(item[0] for item in itertools.takewhile(lambda x: x[1] > 1, c.most_common()))

    def rename_internal_duplicates(node, duplicates):
        for n in traverse(node):
            if n["id"] in duplicates and "children" in n:
                n["id"] = "{}Rep".format(n["id"])

    taxdict = tax.as_dict()
    dupl = get_duplicates(taxdict)
    rename_internal_duplicates(taxdict, dupl)
    res = rec(taxdict) + ";"
    return res


def filter_tree(tree: "ete4.PhyloTree", root_node=None, taxfilter=None):
    newtree = copy.deepcopy(tree)
    if root_node:
        for n in newtree.traverse("preorder"):
            if str(n.name) == str(root_node):
                newtree = n
                break
    if taxfilter:
        for n in newtree.traverse("preorder"):
            if n.name in taxfilter:
                n.detach()
    return newtree


def taxa_index(tree: "ete4.PhyloTree"):
    """build a dictionary from tree-node name to idx"""
    # levelorder is the default traversal of ete3 and ete4; it defines the column order of the profiles
    index = {n.name: i for i, n in enumerate(tree.traverse("levelorder"))}
    return index


def leaf_index(tree: "ete4.PhyloTree"):
    index = {n.name: i for i, n in enumerate(tree.leaves())}
    return index
