import datasketch
import numpy
import pickle
import tables

from .species_tree import clade_ranges, load_species_tree, topmost_nodes_matching
import logging

logger = logging.getLogger(__name__)


def load_forest(root):
    node = root.min_hash_lsh_forest
    if isinstance(node.atom, tables.ObjectAtom):
        # Old format (VLArray of pickled objects)
        return node[0]
    elif isinstance(node.atom, tables.UInt8Atom):
        # New format (EArray of uint8 → bytes → pickle)
        buffer = node[:].tobytes()
        return pickle.loads(buffer)
    else:
        raise TypeError(f"Unsupported atom type: {type(node.atom)}")


class Profiler(object):
    def __init__(self, db):
        root_level_data = db.get_hdf5_handle().get_node("/HOGProfile/ALL")
        self.db = db  # backref
        self.forest = load_forest(root_level_data)
        self.hashes = root_level_data.hashes
        self.species_profile = root_level_data.species_profile
        self.species_tree = load_species_tree(root_level_data)
        self.num_perm = db.get_hdf5_handle().get_node_attr("/HOGProfile/ALL", "num_perm")
        self.species_names = self._build_tree_range_lookup()
        self.tax_range_index = self._select_key_levels_and_ranges()

    def _map_taxid_to_name(self, taxid):
        if int(taxid) == 0:
            return "LUCA"
        tax_row = self.db.tax._taxon_from_numeric(int(taxid))
        return tax_row["Name"].decode()

    def _build_tree_range_lookup(self):
        rows = self.species_tree
        leaves = numpy.flatnonzero(rows["LeafIdx"] >= 0)
        species_names = [""] * len(leaves)
        for i in leaves:
            species_names[rows["LeafIdx"][i]] = self._map_taxid_to_name(rows["TaxonId"][i])
        self._clade_start, self._clade_stop, self._clade_size = clade_ranges(rows)
        return species_names

    def _select_key_levels_and_ranges(self):
        rows = self.species_tree

        def clades(nodes):
            return {
                self._map_taxid_to_name(rows["TaxonId"][i]): (int(self._clade_start[i]), int(self._clade_stop[i]))
                for i in nodes
            }

        indexes = {}
        for size in (50, 200):
            mask = (self._clade_size <= size) & (self._clade_size > 0.1 * size)
            indexes[size] = clades(topmost_nodes_matching(rows, mask))
        indexes["root"] = clades(numpy.flatnonzero(rows["ParentIdx"] == 0))
        return indexes

    def query(self, fam_nr, k=50):
        """Query the database for similar hogs than fam_nr.

        :param int fam_nr: numeric root family
        :param int k: max number of similar root-level hogs being returned"""
        hashvalues = self.hashes[fam_nr].reshape(self.num_perm, 2)
        minhash = datasketch.WeightedMinHash(seed=1, hashvalues=hashvalues)
        similar = self.forest.query(minhash, k=k)

        jaccard_distance = {}
        for sim in similar:
            sval = self.hashes[int(sim)].reshape(self.num_perm, 2)
            shash = datasketch.WeightedMinHash(seed=1, hashvalues=sval)
            jaccard_distance[sim] = shash.jaccard(minhash)

        return ProfileSearchResult(self, fam_nr, similar, jaccard_distance)


class ProfileSearchResult(object):
    def __init__(self, p: Profiler, query_fam, similar_fams, jaccard_distance):
        self.query_fam = query_fam
        self.similar = {int(fam): p.species_profile[int(fam)] for fam in similar_fams}
        self.jaccard_distance = jaccard_distance
        self.query_profile = p.species_profile[query_fam]
        self.tax_classes = p.tax_range_index
        self.species_names = p.species_names
