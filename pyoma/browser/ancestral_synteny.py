import itertools

import networkx as nx
import logging
import numpy
import pyham
import tables
from .tablefmt import AncestralSyntenyRels
from . import db
from .models import Genome
from concurrent.futures import TimeoutError
import multiprocessing as mp
from pebble import ProcessExpired, concurrent

logger = logging.getLogger(__name__)
mp.set_start_method("fork")


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
        try:
            genome = tree_node.genome
        except AttributeError:
            logger.warning("No genome stored for {}".format(tree_node.name))
            continue

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
                # we remove nodes that do not have a parent hog and connect their neighbors
                # (given they have exactly two neighbors) in a working copy that is used
                # to propagate the synteny graph to the higher level.
                G = child.synteny.copy()
                remove_nodes_on_path_without_parents(G)
                # cleanup graph if they get too big
                delete_irrelevant_edges(G, min_importance=7, max_neighbors=8)
                for u, v, weight in G.edges.data("weight", default=1):
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


def remove_nodes_on_path_without_parents(G: nx.Graph):
    """Modifies inplace the Graph such that all nodes
    that have exactly two neighbors and do not have a parent
    HOG are removed and their adjacent nodes are reconnected.

    This works transitively, eg. if A - b - c - D is a graph,
    and b, c have no parents, the pruned graph will be A - D.

    :Example:
    >>> class N:
    ...    def __init__(self, parent=None):
    ...        self.parent = parent
    ...    def __repr__(self):
    ...        return "N({})".format(self.parent)
    >>> A = N('A')
    >>> b = N()
    >>> c = N()
    >>> D = N('D')
    >>> E = N('E')
    >>> G = nx.Graph()
    >>> G.add_weighted_edges_from([(A, b, 1), (b, c, 1), (c, D, 2), (D, E, 1)])
    >>> remove_nodes_on_path_without_parents(G)
    >>> G.edges
    EdgeView([(N(A), N(D)), (N(D), N(E))])

    :param nx.Graph G: input Graph to be analysed
    """
    cnt_removed_nodes = 0
    for node in list(G.nodes):  # convert to list to allow graph modification
        if node.parent is None and len(G.adj[node]) == 2:
            left, right = tuple(G.adj[node])
            if not G.has_edge(left, right):
                G.add_edge(
                    left,
                    right,
                    weight=max(G[left][node]["weight"], G[right][node]["weight"]),
                )
                G.remove_node(node)
                cnt_removed_nodes += 1
    logger.info(
        "removed {} nodes that have synteny but no parent hog.".format(
            cnt_removed_nodes
        )
    )


def delete_irrelevant_edges(G: nx.Graph, min_importance=10, max_neighbors=10):
    """removes edges that are no longer relevant.

    :Example:
    >>> G = nx.Graph()
    >>> G.add_weighted_edges_from([(0, 1, 4), (1, 9, 5)])
    >>> G.add_weighted_edges_from([(1, x, 1) for x in range(2,8)])
    >>> S = G.copy()
    >>> delete_irrelevant_edges(G, min_importance=3)
    >>> G.edges.data('weight')
    EdgeDataView([(0, 1, 4), (1, 9, 5)])
    >>> delete_irrelevant_edges(S, max_neighbors=5, min_importance=5)
    >>> S.edges.data('remove', 0)
    EdgeDataView([(0, 1, 0), (1, 9, 0), (1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1)])

    """
    remove_edge_under_weight = max(min_importance - 2, 1.001)
    cnt_removed_edges = 0
    for u, nbrs in G.adjacency():
        weights = numpy.array([nbr["weight"] for nbr in nbrs.values()])
        weights.sort()
        pos = weights.searchsorted(min_importance)
        to_remove = []
        if pos <= len(weights) - 2:
            # we have at least two well supported edges. remove the ones
            # that are weakly supported.
            for v, eattr in nbrs.items():
                if eattr["weight"] < remove_edge_under_weight:
                    to_remove.append((u, v))
        elif len(weights) > max_neighbors:
            # mark all edge as potential removal, if marked from both sides, remove edge
            for v, eattr in nbrs.items():
                if eattr["weight"] < remove_edge_under_weight:
                    if "remove" in eattr:
                        to_remove.append((u, v))
                    else:
                        eattr["remove"] = 1
        G.remove_edges_from(to_remove)
        cnt_removed_edges += len(to_remove)
    logger.info("removed {} irrelevant edges".format(cnt_removed_edges))


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
    logger.info("starting remove_fork_from_gene_losses...")
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
    logger.info("finished remove_fork_from_gene_losses")


def extract_hog_row_links(graph, hogid_2_hogrow_lookup, level):
    def map_id(id):
        try:
            k = id.index("_")
            query = id[:k]
        except ValueError:
            query = id
        return hogid_2_hogrow_lookup[query]

    edges = []
    nr_errors = 0
    for u, v, w in graph.edges.data("weight", default=1):
        try:
            edges.append((map_id(u), map_id(v), w))
        except KeyError as e:
            nr_errors += 1
            logger.error(
                "unmappable relation ({}, {}) on level {}: {}".format(u, v, level, e)
            )
    if nr_errors > 0:
        logger.warning("{} errors on level {}".format(nr_errors, level))
    return numpy.array(edges, dtype=tables.dtype_from_descr(AncestralSyntenyRels()))


def taxid_from_level(h5, level):
    try:
        taxnode = h5.tax.get_taxnode_from_name_or_taxid(level)
        taxid = int(taxnode["NCBITaxonId"])
    except KeyError as e:
        if level == "LUCA":
            taxid = 0
        else:
            raise
    logger.info("level '{}' -> taxid = {}".format(level, taxid))
    return taxid


@concurrent.process(timeout=1200)
def hogid_2_rownr(h5, level):
    lookup = {}
    level_query = "Level == {!r}".format(level.encode("utf-8"))
    row_iter = h5.get_hdf5_handle().root.HogLevel.where(level_query)
    for row in row_iter:
        lookup[row["ID"].decode()] = row.nrow
    logger.info("hogmap: found {} hog at level {}".format(len(lookup), level))
    return lookup


def level_only_graph_with_hogid_nodes(graph: nx.Graph) -> nx.Graph:
    G = nx.relabel_nodes(graph, lambda n: n.hog_id)
    return G


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
    for cnt, (level, graph) in enumerate(assign_ancestral_synteny(ham), start=1):
        shallow_graph = level_only_graph_with_hogid_nodes(graph)
        future = hogid_2_rownr(h5, level)
        try:
            hogid_lookup = future.result()
            taxid = taxid_from_level(h5, level)
            edges = extract_hog_row_links(shallow_graph, hogid_lookup, level)
            synteny_graphs_mapped[taxid] = edges
        except TimeoutError as error:
            logger.error(
                "hogid_2_rownr took longer than %d seconds".format(error.args[1])
            )
        except ProcessExpired as error:
            logger.error("{}. Exit code: {}}".format(error, error.exitcode))
        except Exception as error:
            logger.exception("hogid_2_rownr raised {}".format(error))

        if cnt % 200 == 0:
            logger.info("closing and reopening {} file".format(h5name))
            h5.close()
            h5 = db.Database(h5name)
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
