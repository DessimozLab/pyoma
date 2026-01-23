import logging
import pickle
import collections

import tables
import pandas
import rapidfuzz
from omataxonomy import Taxonomy

from ...common import auto_open

logger = logging.getLogger(__name__)


def _load_tax_mappings(tax):
    gtdb2ncbi = collections.defaultdict(list)
    c = tax.db.execute('SELECT * from synonym where spname like "ncbi%";')
    for gtdb, ncbi_text in c.fetchall():
        ncbi = int(ncbi_text[ncbi_text.index(":") + 1 :])
        gtdb2ncbi[gtdb].append(ncbi)

    ncbi2gtdb = collections.defaultdict(list)
    for g, n in gtdb2ncbi.items():
        for k in n:
            ncbi2gtdb[k].append(g)

    return gtdb2ncbi, ncbi2gtdb


def _load_gs_data(h5_path: str) -> pandas.DataFrame:
    with tables.open_file(h5_path) as db:
        gs = pandas.DataFrame(db.get_node("/Genome").read())
        tax = pandas.DataFrame(db.get_node("/Taxonomy").read())
    gs["GenomeId"] = gs["NCBITaxonId"]
    special = (gs["GenomeId"] < 0) & (gs["GenomeId"] > -(len(gs) + 100))
    gs = pandas.merge(gs, tax, left_on="NCBITaxonId", right_on="NCBITaxonId")
    gs.loc[special, "NCBITaxonId"] = gs.loc[special, "ParentTaxonId"]
    for col in gs.columns:
        if gs[col].dtype == "object":
            gs[col] = gs[col].str.decode("utf-8")
    return gs


def extract_source(row: pandas.Series) -> str:
    if "ensembl.org" in row["Source"]:
        return "Ensembl"
    elif "ensemblgenomes.org" in row["Source"]:
        return "EnsemblGenomes"
    elif "/GCF/" in row["Source"]:
        return "RefSeq"
    elif "/GCA/" in row["Source"]:
        return "Genbank"
    elif "ncbi.nlm.nih.gov" in row["Source"]:
        return "NCBI"
    elif "genome_reviews" in row["Source"]:
        return "GenomeReviews"
    elif "ebi.ac.uk" in row["Source"]:
        return "EBI"
    elif "flybase.net" in row["Source"]:
        return "FlyBase"
    elif "jgi.doe.gov" in row["Source"] or "jgi-psf.org" in row["Source"]:
        return "JGI"
    elif "wormbase.org" in row["Url"]:
        return "WormBase"
    elif "silkdb.org" in row["Url"]:
        return "SilkDB"
    elif "dictybase.org" in row["Url"]:
        return "DictyBase"
    elif "phytozome" in row["Url"] or "phytozome" in row["Source"]:
        return "Phytozome"
    elif "broad.mit.edu" in row["Url"] or "broad" in row["Source"].lower():
        return "BROAD Institute"
    elif "modENCODE" in row["Release"] or "modENCODE" in row["Source"]:
        return "modENCODE"
    elif "ugent.be/plaza" in row["Url"]:
        return "Plaza"
    elif "ugent.be/orcae" in row["Url"]:
        return "OrcAE"
    elif "cottongen.org" in row["Url"]:
        return "CottonGen"
    elif "ginkgo." in row["Url"]:
        return "GinkgoDB"
    else:
        parts = row["Release"].split("; ")
        if len(parts) > 1:
            return parts[0]
        return "n/a"


def add_gtdb_and_ncbi_columns(gs: pandas.DataFrame, tax: Taxonomy, gtdb2ncbi: dict):
    taxids = list(gs["NCBITaxonId"])
    tax2acc = dict(zip(taxids, tax.translate_to_names(taxids)))
    gs = gs.copy(deep=True)

    gs["GTDB_Accession"] = "n/a"
    gs["GTDB_Accession"] = gs["GTDB_Accession"].where(gs["NCBITaxonId"] >= 0, gs["NCBITaxonId"].map(tax2acc))
    gs["NCBITaxonId"] = gs["NCBITaxonId"].where(
        gs["NCBITaxonId"] >= 0, gs["NCBITaxonId"].map(lambda x: gtdb2ncbi.get(x, [0])[0])
    )
    return gs


def map_metadata_columns(gs: pandas.DataFrame) -> pandas.DataFrame:
    gs = gs.copy(deep=True)
    if {"Source", "Release", "Url"} - set(gs.columns):
        logger.info("not all metadata columns present to map Source")
        gs["SourceInfo"] = gs["Source"] if "Source" in gs.columns else "n/a"
        if "Release" not in gs.columns:
            gs["Release"] = "n/a"
    else:
        gs["SourceInfo"] = gs.apply(extract_source, axis=1)
    return gs


def create_treenode_to_genome_map(tree, tax2code):
    res = {}
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            try:
                spp = set(tax2code[int(n.name)])
            except KeyError:
                spp = set()
            n.add_feature("spset", spp)
            res[int(n.name)] = n.spset
        else:
            spset = set()
            for x in n.get_children():
                spset |= x.spset
            n.add_feature("spset", spset)
            res[int(n.name)] = spset
    return res


def find_top_overlaping_nodes(source, target):
    max_pairs = {}
    for src_id, src_set in source.items():
        res = []
        for tgt_id, tgt_set in target.items():
            jac = (
                2 * len(tgt_set.intersection(src_set)) / (len(tgt_set) + len(src_set))
                if (len(tgt_set) + len(src_set)) > 0
                else 0.0
            )
            res.append((tgt_id, jac))
        res.sort(key=lambda x: -x[1])
        max_pairs[src_id] = res[:50]
    return max_pairs


def _filter_candidates(tax, cand, name, threshold=0.4):
    for k in range(len(cand)):
        if cand[k][1] < threshold * cand[0][1]:
            break
    rel_cand = [
        (
            x[0],
            tax.translate_to_names([x[0]])[0],
            x[1],
            rapidfuzz.fuzz.ratio(tax.translate_to_names([x[0]])[0], name),
        )
        for x in cand[:k]
    ]
    rel_cand.sort(key=lambda e: (-e[2], -e[3]))
    return rel_cand


def collect_gtdb_map_info(cand_map, tax2code, code2tax, tax, pyomatax, tree):
    ancestral, extant = [], []
    for gtdb_id, cand in cand_map.items():
        if gtdb_id in tax2code:
            # gtdb id is extant species
            this_one = {"gtdb_taxid": gtdb_id, "gtdb_acc": tax.translate_to_names([gtdb_id])[0]}
            ncbi_this = [xx for c in tax2code[gtdb_id] for xx in code2tax[c] if xx > 0]
            if len(ncbi_this) > 0:
                ncbi_this = ncbi_this[0]
                this_one["ncbi_taxid"] = ncbi_this
                this_one["sciname"] = tax.translate_to_names([ncbi_this])[0]
            extant.append(this_one)
        else:
            this_one = {"gtdb_taxid": gtdb_id, "gtdb_name": tax.translate_to_names([gtdb_id])[0]}
            if this_one["gtdb_name"] in ("d__Bacteria", "d__Archaea"):
                this_one["valid_oma_tax"] = False
                if this_one["gtdb_name"] == "d__Bacteria":
                    this_one["representative"] = 2
                elif this_one["gtdb_name"] == "d__Archaea":
                    this_one["representative"] = 2157
                ancestral.append(this_one)
                continue

            if gtdb_id not in pyomatax["NCBITaxonId"]:
                this_one["valid_oma_tax"] = False
                tnode = tree.search_nodes(taxid=gtdb_id)[0]
                while tnode.taxid not in pyomatax["NCBITaxonId"]:
                    tnode = tnode.children[0]
                this_one["representative"] = tnode.taxid
                logger.info(
                    f"{gtdb_id} - {this_one['gtdb_name']} not in pyomatax -> {tnode.taxid} - {tax.translate_to_names([tnode.taxid])[0]}"
                )
            else:
                this_one["valid_oma_tax"] = True
            this_one["ncbi_candidates"] = _filter_candidates(tax, cand, this_one["gtdb_name"])
            ancestral.append(this_one)
    return ancestral, extant


def collect_ncbi_map_info(cand_map, tax2code, code2tax, tax, pyomatax, tree):
    ancestral, extant = [], []
    for ncbi_id, cand in cand_map.items():
        if ncbi_id in tax2code:
            # ncbi id is extant species
            this_one = {"ncbi_taxid": ncbi_id, "sciname": tax.translate_to_names([ncbi_id])[0]}
            gtdb_this = [xx for code in tax2code[ncbi_id] for xx in code2tax[code] if xx < 0]
            if len(gtdb_this) < 2:
                if len(gtdb_this) == 1:
                    gtdb_id = gtdb_this[0]
                    this_one["gtdb_taxid"] = gtdb_id
                    this_one["gtdb_acc"] = tax.translate_to_names([gtdb_id])[0]
                extant.append(this_one)
                continue

        # either not a leaf or a leaf that maps to multiple gtdb taxa (--> same as internal node)
        this_one = {"ncbi_taxid": ncbi_id, "sciname": tax.translate_to_names([ncbi_id])[0]}

        if 2759 in tax.get_lineage(ncbi_id) or ncbi_id in (2, 2157):
            if ncbi_id not in pyomatax["NCBITaxonId"]:
                this_one["valid_oma_tax"] = False
                tnode = tree.search_nodes(taxid=ncbi_id)[0]
                while tnode.taxid not in pyomatax["NCBITaxonId"]:
                    try:
                        tnode = tnode.children[0]
                    except IndexError:
                        logger.info(f"Cannot find a genome for {ncbi_id} - {this_one['sciname']}")
                        break
                if tnode.taxid not in pyomatax["NCBITaxonId"]:
                    continue
                this_one["representative"] = tnode.taxid
                logger.info(
                    f"{ncbi_id} - {this_one['sciname']} not in pyomatax -> {tnode.taxid} - {tax.translate_to_names([tnode.taxid])[0]}"
                )
            else:
                this_one["valid_oma_tax"] = True
            ancestral.append(this_one)
            continue  # eukaryotic genomes do not have a gtdb mapping
        this_one["gtdb_candidates"] = _filter_candidates(tax, cand, this_one["sciname"])
        ancestral.append(this_one)
    return ancestral, extant


def iter_extant_unique(lst):
    seen = set()
    for e in lst:
        t = (e.get("ncbi_taxid"), e.get("gtdb_taxid"))
        if t not in seen:
            seen.add(t)
            yield e


def map_gtdb_and_ncbi(tax: Taxonomy, gs: pandas.DataFrame, db_path: str, gtdb2ncbi: dict) -> dict:
    gtdb_mask = gs["NCBITaxonId"] < 0
    oma_gtdb = list(gs.loc[gtdb_mask, "NCBITaxonId"])
    oma_ncbi = list(gtdb2ncbi.get(id_, [0])[0] for id_ in oma_gtdb)
    oma_ncbi.extend(list(gs.loc[~gtdb_mask, "NCBITaxonId"]))
    tax2code = collections.defaultdict(set)
    for row in gs.itertuples(index=False):
        tax2code[row.NCBITaxonId].add(row.UniProtSpeciesCode)
        if row.NCBITaxonId < 0:
            tax2code[gtdb2ncbi.get(row.NCBITaxonId, 0)[0]].add(row.UniProtSpeciesCode)

    tree_gtdb = tax.get_topology(oma_gtdb, intermediate_nodes=True, annotate=True)
    tree_ncbi = tax.get_topology(oma_ncbi, intermediate_nodes=True, annotate=True)

    gtdb_nodes_to_genome_set = create_treenode_to_genome_map(tree_gtdb, tax2code)
    ncbi_nodes_to_genome_set = create_treenode_to_genome_map(tree_ncbi, tax2code)
    gtdb_to_ncbi_cand = find_top_overlaping_nodes(gtdb_nodes_to_genome_set, ncbi_nodes_to_genome_set)
    ncbi_to_gtdb_cand = find_top_overlaping_nodes(ncbi_nodes_to_genome_set, gtdb_nodes_to_genome_set)

    code2tax = collections.defaultdict(set)
    for k, codes in tax2code.items():
        for code in codes:
            code2tax[code].add(k)

    with tables.open_file(db_path, "r") as h5:
        taxtab = h5.get_node("/Taxonomy").read()
    ancestral, extant = collect_gtdb_map_info(gtdb_to_ncbi_cand, tax2code, code2tax, tax, taxtab, tree_gtdb)
    ncbi_maps = collect_ncbi_map_info(ncbi_to_gtdb_cand, tax2code, code2tax, tax, taxtab, tree_ncbi)
    ancestral.extend(ncbi_maps[0])
    extant.extend(ncbi_maps[1])

    map_data = {"ancestral": ancestral, "extant": list(iter_extant_unique(extant))}
    return map_data


def produce_mapping_and_download_files(db_path, tax_path, download_path, mapping_path):
    taxonomy = Taxonomy(tax_path)
    gtdb2ncbi, ncbi2gtdb = _load_tax_mappings(taxonomy)
    gs_basic = _load_gs_data(db_path)
    gs_mapped = add_gtdb_and_ncbi_columns(gs_basic, taxonomy, gtdb2ncbi)
    gs_mapped = map_metadata_columns(gs_mapped)
    gs_mapped = gs_mapped[
        ["UniProtSpeciesCode", "GenomeId", "NCBITaxonId", "SourceInfo", "SciName", "GTDB_Accession", "Release"]
    ]
    gs_mapped.rename(
        columns={
            "UniProtSpeciesCode": "OMA_Code",
            "GenomeId": "OMA_Taxon_ID",
            "NCBITaxonId": "NCBI_Taxon_ID",
            "SourceInfo": "Source",
            "SciName": "Scientific_Name",
        },
        inplace=True,
    )
    with auto_open(download_path, "wt") as fh:
        fh.write("# Mapping of OMA species codes to NCBI taxon IDs and scientific names.\n")
        fh.write("# Note: OMA species codes are whenever possible identical to UniProt codes.\n")
        fh.write(
            "# As OMA Taxon IDs, the NCBI Taxon IDs are used for eukaryotic genomes, and negative\n"
            "# values based on the GTDB Accessions for prokaryotic genomes.\n"
        )
        fh.write(
            "# Format: OMA_Code<tab>OMA_Taxon_ID<tab>NCBI_Taxon_ID<tab>GTDB_Accession<tab>Scientific_Name<tab>Source<tab>Release/Version\n"
        )
        gs_mapped[
            ["OMA_Code", "OMA_Taxon_ID", "NCBI_Taxon_ID", "GTDB_Accession", "Scientific_Name", "Source", "Release"]
        ].sort_values("OMA_Code").to_csv(fh, sep="\t", index=False, na_rep="n/a")

    mapping_data = map_gtdb_and_ncbi(taxonomy, gs_basic, db_path, gtdb2ncbi)
    with auto_open(mapping_path, "wb") as fh:
        pickle.dump(mapping_data, fh)
