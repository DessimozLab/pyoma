"""Storage of the species tree used by the HOG profiles.

At query time only the layout of the profile columns is needed, i.e. which
species a column of ``species_profile`` corresponds to and which (contiguous)
column ranges the clades span. The tree is therefore stored as a table at
``/HOGProfile/ALL/species_tree`` with one row per node in preorder (see
:class:`pyoma.browser.tablefmt.ProfileSpeciesTreeTable`). This keeps the
hdf5 file independent of any tree library.

Older builds stored a pickled ete3 tree (VLArray with ObjectAtom). These are
converted to the table format on load without requiring ete3.
"""

import io
import pickle
from typing import Any, Callable, Iterable, Tuple

import numpy
import tables

from ..tablefmt import ProfileSpeciesTreeTable

species_tree_dtype = tables.dtype_from_descr(ProfileSpeciesTreeTable)


def taxid_from_node_name(name) -> int:
    """converts a node name of the profile species tree to a taxid.

    Internal nodes whose taxid also appears as leaf carry a 'Rep' suffix."""
    name = str(name)
    if name.endswith("Rep"):
        name = name[: -len("Rep")]
    # at some point the taxid's were uint32 with negative values
    # need to cast them to int32
    return (int(name) + 2**31) % 2**32 - 2**31


def tree_to_rows(root, children_of: Callable[[Any], Iterable], name_of: Callable[[Any], Any]) -> numpy.ndarray:
    """converts a tree into the rows of the species tree table.

    Nodes are listed in preorder and leaves are numbered in the same order,
    which corresponds to the order of the columns of the species profiles.

    :param root: root node of the tree
    :param children_of: function returning the children of a node
    :param name_of: function returning the name of a node"""
    rows = []
    nr_leaves = 0
    stack = [(root, -1)]
    while stack:
        node, parent = stack.pop()
        children = list(children_of(node))
        if children:
            leaf_idx = -1
        else:
            leaf_idx, nr_leaves = nr_leaves, nr_leaves + 1
        rows.append((taxid_from_node_name(name_of(node)), parent, leaf_idx))
        idx = len(rows) - 1
        stack.extend((child, idx) for child in reversed(children))
    return numpy.array(rows, dtype=species_tree_dtype)


def ete_tree_to_rows(tree) -> numpy.ndarray:
    """converts an ete tree (ete3 or ete4) into the species tree table rows"""
    return tree_to_rows(tree, children_of=lambda n: n.children, name_of=lambda n: n.name)


def store_species_tree(h5: tables.File, where, rows: numpy.ndarray, **kwargs) -> tables.Table:
    return h5.create_table(where, "species_tree", obj=rows, expectedrows=len(rows), **kwargs)


class _LegacyEteNode(object):
    """placeholder for the ete classes referenced in legacy pickles"""


class _LegacyEteUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module.split(".")[0] in ("ete2", "ete3"):
            return _LegacyEteNode
        raise pickle.UnpicklingError("unexpected global '{}.{}' in pickled species tree".format(module, name))


def load_legacy_pickled_tree(data: bytes) -> numpy.ndarray:
    """converts a pickled ete3 tree into the species tree table rows without requiring ete3."""
    legacy_root = _LegacyEteUnpickler(io.BytesIO(data)).load()
    return tree_to_rows(
        legacy_root,
        children_of=lambda n: vars(n).get("_children", []),
        name_of=lambda n: vars(n).get("name", vars(n).get("_name")),
    )


def load_species_tree(root) -> numpy.ndarray:
    """loads the species tree table rows from the HOGProfile group ``root``."""
    node = root.species_tree
    if isinstance(node, tables.Table):
        return node.read()
    elif isinstance(node, tables.VLArray) and isinstance(node.atom, tables.ObjectAtom):
        # Old format (pickled ete3 tree). Read the raw bytes to avoid
        # that pytables unpickles it (which would require ete3).
        return load_legacy_pickled_tree(node._read_array(0, 1, 1)[0].tobytes())
    else:
        raise TypeError(f"Unsupported species_tree node: {node}")


def clade_ranges(rows: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """computes for every node the range of profile columns its clade covers.

    :returns: arrays start, stop (exclusive) and size (number of leaves)"""
    is_leaf = rows["LeafIdx"] >= 0
    start = numpy.where(is_leaf, rows["LeafIdx"], numpy.iinfo(numpy.int32).max).astype(numpy.int64)
    stop = numpy.where(is_leaf, rows["LeafIdx"] + 1, 0).astype(numpy.int64)
    size = is_leaf.astype(numpy.int64)
    # in preorder, children come after their parent, so a reverse scan
    # has processed all children before reaching their parent.
    for i in range(len(rows) - 1, 0, -1):
        p = rows["ParentIdx"][i]
        start[p] = min(start[p], start[i])
        stop[p] = max(stop[p], stop[i])
        size[p] += size[i]
    return start, stop, size


def topmost_nodes_matching(rows: numpy.ndarray, mask: numpy.ndarray) -> numpy.ndarray:
    """returns the indices (in preorder) of the nodes in mask that have no
    ancestor in mask, i.e. a traversal from the root stops at these nodes."""
    covered = numpy.zeros(len(rows), dtype=bool)
    selected = []
    for i in range(len(rows)):
        p = rows["ParentIdx"][i]
        if p >= 0 and covered[p]:
            covered[i] = True
        elif mask[i]:
            covered[i] = True
            selected.append(i)
    return numpy.array(selected, dtype=numpy.int64)
