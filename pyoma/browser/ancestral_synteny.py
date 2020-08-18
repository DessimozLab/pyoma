import networkx as nx
import logging
import pyham
from . import db
from .models import Genome

logger = logging.getLogger(__name__)


def assign_extant_syteny(h5, ham):
    genomes = [Genome(h5, g) for g in h5.get_hdf5_handle().root.Genome.read()]
    for g in genomes:
        proteins = h5.main_isoforms(g.uniprot_species_code)
        proteins.sort(order=["Chromosome", "LocusStart"])
        graph = nx.Graph()
        gene1 = ham.get_gene_by_id(proteins[0]["EntryNr"])
        graph.add_node(gene1)
        for i in range(1, len(proteins)):
            gene = ham.get_gene_by_id(proteins[i]["EntryNr"])
            graph.add_node(gene)
            if proteins[i - 1]["Chromosome"] == proteins[i]["Chromosome"]:
                # connect genes only if they are on the same chromosome
                graph.add_edge(gene1, gene, weight=1)
            gene1 = gene
        if logger.isEnabledFor(logging.INFO):
            logger.info(
                "Synteny for extant genome {} established. |V|={}, |E|={}, |CC|={}".format(
                    g.uniprot_species_code,
                    graph.number_of_nodes(),
                    graph.number_of_edges(),
                    nx.number_connected_components(graph),
                )
            )
        gene.genome.taxon.add_feature("synteny", graph)


def assign_ancestral_synteny(ham):
    for tree_node in ham.taxonomy.tree.traverse("postorder"):
        if isinstance(tree_node.genome, pyham.ExtantGenome):
            logger.info(
                "synteny for extant genome '{}' already assigned".format(
                    tree_node.genome
                )
            )
            continue
        elif isinstance(tree_node.genome, pyham.AncestralGenome):
            graph = nx.Graph()
            graph.add_nodes_from(tree_node.genome.genes)
            for child in tree_node.children:
                for u, v, weight in child.synteny.edges.data("weight", default=1):
                    if u.parent is not None and v.parent is not None:
                        assert u.parent in graph.nodes
                        assert v.parent in graph.nodes
                        parent_edge = (u.parent, v.parent)
                        if graph.has_edge(*parent_edge):
                            graph[parent_edge[0]][parent_edge[1]]["weight"] += weight
                        else:
                            graph.add_edge(*parent_edge, weight=weight)
            tree_node.add_feature("synteny", graph)
            logger.info(
                "Synteny for ancestral genome '{}' created. |V|={}, |E|={}, |CC|={}".format(
                    tree_node.name,
                    graph.number_of_nodes(),
                    graph.number_of_edges(),
                    nx.number_connected_components(graph),
                )
            )


def infer_synteny(orthoxml, h5name, tree):
    h5 = db.Database(h5name)
    ham = pyham.Ham(
        tree_file=tree, hog_file=orthoxml, tree_format="newick", use_internal_name=True
    )
    assign_extant_syteny(h5, ham)
    assign_ancestral_synteny(ham)


if __name__ == "__main__":
    orthoxml = "/Volumes/TOSHIBA EXT/basf/All.May2020/downloads/oma-hogs.orthoXML.gz"
    tree = "/Volumes/TOSHIBA EXT/basf/All.May2020/downloads/speciestree.nwk"
    h5db = "/Volumes/TOSHIBA EXT/basf/All.May2020/data/OmaServer.h5"
    infer_synteny(orthoxml, h5db, tree)
