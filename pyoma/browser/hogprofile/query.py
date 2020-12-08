import datasketch
import numpy

from .tree_helper import leaf_index
import logging

logger = logging.getLogger(__name__)


class Profiler(object):
    def __init__(self, db):
        root_level_data = db.get_hdf5_handle().get_node("/HOGProfile/ALL")
        self.db = db  # backref
        self.forest = root_level_data.min_hash_lsh_forest[0]
        self.hashes = root_level_data.hashes
        self.species_profile = root_level_data.species_profile
        self.species_tree = root_level_data.species_tree[0]
        self.num_perm = db.get_hdf5_handle().get_node_attr(
            "/HOGProfile/ALL", "num_perm"
        )
        self.species_names = self._build_tree_range_lookup()
        self.tax_range_index = self._select_key_levels_and_ranges()

    def _map_taxid_to_name(self, taxid):
        if int(taxid) == 0:
            return "LUCA"
        tax_row = self.db.tax._taxon_from_numeric(int(taxid))
        return tax_row["Name"].decode()

    def _build_tree_range_lookup(self):
        leaf_idx = leaf_index(self.species_tree)
        species_names = [""] * len(leaf_idx)
        for x in self.species_tree.iter_leaves():
            species_names[leaf_idx[x.name]] = self._map_taxid_to_name(x.name)

        def idx_map(node):
            try:
                taxid = int(node.name)
            except ValueError:
                taxid = int(node.name[: node.name.index("Rep")])
            # at some point the taxid's were uint32 with negative values
            # need to cast them via np.int32
            taxid = int(numpy.int32(taxid))
            sciname = self._map_taxid_to_name(taxid)
            if node.is_leaf():
                idx = leaf_idx[node.name]
                node.add_features(range=(idx, idx + 1), size=1, sciname=sciname)
            else:
                min_range, max_range, size = 10000, 0, 0
                for child in node.children:
                    idx_map(child)
                    if child.range[0] < min_range:
                        min_range = child.range[0]
                    if child.range[1] > max_range:
                        max_range = child.range[1]
                    size += child.size
                node.add_features(
                    range=(min_range, max_range), size=size, sciname=sciname
                )

        # annotate speciestree with range features
        idx_map(self.species_tree)
        return species_names

    def _select_key_levels_and_ranges(self):
        def nodes_of_size(node, size):
            for n in node.iter_leaves(is_leaf_fn=lambda n: size >= n.size > 0.1 * size):
                yield n

        indexes = {}
        for size in (50, 200):
            indexes[size] = {
                n.sciname: n.range for n in nodes_of_size(self.species_tree, size)
            }
        indexes["root"] = {n.sciname: n.range for n in self.species_tree.children}
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
