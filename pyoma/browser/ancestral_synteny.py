import itertools

import networkx as nx
import logging

import numpy
import pyham
import tables

from .tablefmt import AncestralSyntenyRels
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
        gene1.genome.taxon.add_feature("synteny", graph)


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
                    if (
                        u.parent is not None
                        and v.parent is not None
                        and u.parent != v.parent
                    ):
                        assert u.parent in graph.nodes
                        assert v.parent in graph.nodes
                        parent_edge = (u.parent, v.parent)
                        if graph.has_edge(*parent_edge):
                            graph[parent_edge[0]][parent_edge[1]]["weight"] += weight
                        else:
                            graph.add_edge(*parent_edge, weight=weight)
            logger.info("build graph for {}: {}".format(tree_node.name, nx.info(graph)))
            remove_forks_from_gene_losses(graph)
            tree_node.add_feature("synteny", graph)
            logger.info(
                "Synteny for ancestral genome '{}' created. |V|={}, |E|={}, |CC|={}".format(
                    tree_node.name,
                    graph.number_of_nodes(),
                    graph.number_of_edges(),
                    nx.number_connected_components(graph),
                )
            )
            yield tree_node.name, graph


def remove_forks_from_gene_losses(G: nx.Graph):
    """Remove forks in the graph that originate from gene losses
    on the child level.

    <h_a>  --  <h_i>  --  <h_j>  -- <h_b>
                 \           /
                  -- <h_k> --

    this is the base case to remove the edge between h_i and h_j.
    <h_k> has to be at least one hog, but can also be more, but
    the path from h_i to h_j via h_k needs to be simple, e.g. no
    further branching.
    """
    deg_3_nodes = (n[0] for n in G.degree if n[1] == 3)
    pot_cases = list(
        filter(lambda uv: G.has_edge(*uv), itertools.combinations(deg_3_nodes, 2))
    )
    removed_forks = 0
    for u, v in pot_cases:
        pgen = nx.shortest_simple_paths(G, u, v)
        direct = next(pgen)  # direct edge between u,v. to be removed
        try:
            next_path = next(pgen)
        except StopIteration:
            # no other path, continue with next pot_case
            continue
        try:
            more_path = next(pgen)
            continue  # we don't want cases that have more possible paths
        except StopIteration:
            pass

        # remove direct path, reweight next_path
        w = G.edges[u, v]["weight"]
        G.remove_edge(u, v)
        for i in range(len(next_path) - 1):
            G.edges[next_path[i], next_path[i + 1]]["weight"] += w
        removed_forks += 1
        logger.debug(
            "removed edge from ({},{}), found other path of len {}: {}".format(
                u, v, len(next_path), next_path
            )
        )
    logger.info(
        "removed {} forks in G(|V|={},|E|={})".format(
            removed_forks, G.number_of_nodes(), G.number_of_edges()
        )
    )


def extract_hog_row_links(h5, level, graph):
    lookup = {}
    for row in h5.get_hdf5_handle().root.HogLevel.where(
        "Level == {!r}".format(level.encode("utf-8"))
    ):
        lookup[row["ID"].decode()] = row.nrow

    try:
        taxnode = h5.tax.get_taxnode_from_name_or_taxid(level)
        taxid = int(taxnode["NCBITaxonId"])
    except KeyError as e:
        if level == "LUCA":
            taxid = 0
        else:
            raise

    def map_id(id):
        try:
            k = id.index("_")
            query = id[:k]
        except ValueError:
            query = id
        return lookup[query]

    edges = []
    nr_errors = 0
    for u, v, w in graph.edges.data("weight", default=1):
        try:
            edges.append((map_id(u.hog_id), map_id(v.hog_id), w))
        except KeyError as e:
            nr_errors += 1
            logger.error(
                "unmappable relation ({}, {}) on level {}: {}".format(u, v, level, e)
            )
    if nr_errors > 0:
        logger.warning("{} errors on level {}".format(nr_errors, level))
    return (
        taxid,
        numpy.array(edges, dtype=tables.dtype_from_descr(AncestralSyntenyRels())),
    )


def infer_synteny(orthoxml, h5name, tree):
    h5 = db.Database(h5name)
    ham = pyham.Ham(
        tree_file=tree,
        hog_file=orthoxml,
        tree_format="newick",
        use_internal_name=True,
        species_resolve_mode="OMA",
    )
    assign_extant_syteny(h5, ham)
    synteny_graphs_mapped = {}
    for level, graph in assign_ancestral_synteny(ham):
        taxid, edges = extract_hog_row_links(h5, level, graph)
        synteny_graphs_mapped[taxid] = edges
    h5.close()
    return ham, synteny_graphs_mapped


def plot_graph_neighborhood(G, node, steps=2):
    import matplotlib.pyplot as plt

    nodes = [node] + [v for u, v in nx.bfs_edges(G, source=node, depth_limit=steps)]
    S = G.subgraph(nodes)
    colors = [w for u, v, w in S.edges.data("weight", 1)]
    weights = [w for u, v, w in S.edges.data("weight", 1)]
    options = {
        "node_color": "#A0CBE2",
        "edge_color": colors,
        "width": weights,
        "edge_cmap": plt.cm.jet,
        "with_labels": True,
    }
    pos = nx.spring_layout(S, iterations=200)
    nx.draw(S, pos=pos, **options)
    labels = nx.get_edge_attributes(S, "weight")
    nx.draw_networkx_edge_labels(S, pos=pos, edge_labels=labels)
    plt.show()


def write_syntenygraphs_to_hdf5(h5name, hog_edges_per_level):
    with tables.open_file(h5name, "a") as h5:
        node = h5.create_group(
            "/",
            "AncestralSynteny",
            title="Graph of HOGs (identified by rowidx in HogLevel table). Per tax level",
        )
        for level, edges in hog_edges_per_level.items():
            h5.create_table(
                node, "tax{}".format(level), obj=edges, expectedrows=len(edges)
            )


if __name__ == "__main__":
    orthoxml = "/Volumes/TOSHIBA EXT/basf/All.May2020/downloads/oma-hogs.orthoXML.gz"
    tree = "/Volumes/TOSHIBA EXT/basf/All.May2020/downloads/speciestree.nwk"
    h5db = "/Volumes/TOSHIBA EXT/basf/All.May2020/data/OmaServer.h5"
    infer_synteny(orthoxml, h5db, tree)
