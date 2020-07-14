import ete3
import numpy as np
import pandas as pd


def generate_treeweights(tree: ete3.PhyloNode, taxa_index):
    # weighing function for tax level, masking levels etc. sets all weights to 1
    """
    Generate the weights of each taxonomic level to be applied during the
    constructin of weighted minhashes
    :param tree: trimmed species tree
    :param taxa_index: dict mapping taxa to columns
    :return: weights: a vector of weights for each tax level
    """

    weights = np.zeros((3 * len(taxa_index), 1))
    for i in range(3):
        for n in tree.traverse():
            weights[len(taxa_index) * i + taxa_index[n.name]] = 1
    return weights


def hash_tree_profile(tp: ete3.PhyloNode, taxa_index, species_index, treeweights, wmg):
    """
    Generate a weighted minhash and binary matrix row for a tree profile

    :param tp: a pyham tree profile
    :param taxa_index: dict mapping taxa to columns
    :param treeweights: a vector of weights for each tax levels
    :param wmg: Datasketch weighted minhash generator
    :return hog_matrix: a vector of weights for each tax level
    :return weighted_hash: a weighted minhash of a HOG

    """

    losses = [
        taxa_index[n.name] for n in tp.traverse() if n.lost and n.name in taxa_index
    ]
    dupl = [
        taxa_index[n.name] for n in tp.traverse() if n.dupl and n.name in taxa_index
    ]
    presence = [
        taxa_index[n.name]
        for n in tp.traverse()
        if n.nbr_genes > 0 and n.name in taxa_index
    ]
    covered_species = [
        species_index[n.name]
        for n in tp.iter_leaves()
        if n.nbr_genes > 0 and n.name in species_index
    ]
    indices = dict(zip(["presence", "loss", "dup"], [presence, losses, dupl]))
    hog_matrix_weighted = np.zeros((1, 3 * len(taxa_index)))
    hog_matrix_binary = np.zeros((1, 3 * len(taxa_index)))
    species_matrix = np.zeros(len(species_index))
    species_matrix[covered_species] = 1
    for i, event in enumerate(indices):
        if len(indices[event]) > 0:
            hogindex = np.asarray(indices[event]) + i * len(taxa_index)
            hog_matrix_weighted[:, hogindex] = treeweights[hogindex, :].ravel()
            hog_matrix_binary[:, hogindex] = 1
    weighted_hash = wmg.minhash(list(hog_matrix_weighted.flatten()))
    return hog_matrix_binary, weighted_hash, species_matrix


def row2hash(row, taxa_index, species_index, treeweights, wmg):
    """
    turn a dataframe row with an orthoxml file to hash and matrix row
    :param row: lsh builder dataframe row (with family, treeprofile)
    :param taxa_index: dict mapping taxa to columnsfam
    :param treeweights: a vector of weights for each tax levels
    :param wmg: Datasketch weighted minhash generator
    :return: hog_matrix: a vector of weights for each tax level
    :return: weighted_hash: a weighted minhash of a HOG
    """
    # convert a dataframe row to a weighted minhash
    fam, treemap = row.tolist()
    hog_matrix, weighted_hash, species_matrix = hash_tree_profile(
        treemap, taxa_index, species_index, treeweights, wmg
    )
    return pd.Series(
        [weighted_hash, hog_matrix, species_matrix], index=["hash", "rows", "species"]
    )
