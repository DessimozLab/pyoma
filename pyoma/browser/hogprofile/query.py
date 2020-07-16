import datasketch


class Profiler(object):
    def __init__(self, db):
        root_level_data = db.get_hdf5_handle().get_node("/HOGProfile/ALL")
        self.db = db  # backref
        self.forest = root_level_data.min_hash_lsh_forest[0]
        self.hashes = root_level_data.hashes
        self.species_profile = root_level_data.species_profile
        self.num_perm = db.get_hdf5_handle().get_node_attr(
            "/HOGProfile/ALL", "num_perm"
        )

    def query(self, fam_nr, k=50):
        """Query the database for similar hogs than fam_nr.

        :param int fam_nr: numeric root family
        :param int k: max number of similar rootlevel hogs being returned"""
        hashvalues = self.hashes[fam_nr].reshape(self.num_perm, 2)
        minhash = datasketch.WeightedMinHash(seed=1, hashvalues=hashvalues)
        similar = self.forest.query(minhash, k=k)
        similar_fams = {int(fam): self.species_profile[int(fam)] for fam in similar}
        return similar_fams
