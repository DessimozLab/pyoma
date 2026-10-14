import collections
import logging
import os
import resource
import sys
import warnings
from contextlib import ExitStack
from collections import defaultdict
from textwrap import dedent
from typing import Tuple

import numpy as np
import tables
from tables import PerformanceWarning
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from tqdm import tqdm

from .main import setup_logging
from ..db import Database
from ..hoghelper import are_orthologous
from ..models import ProteinEntry
from ..build.taxonomy_alignment import produce_mapping_and_download_files
from ..linkout import ncbi as ncbi_linkout
from ..linkout import uniprot as uniprot_linkout
from ...common import auto_open

logger = logging.getLogger(__name__)


def open_xref_outputs(args):
    """Return dict: SourceXRef -> open file handle"""
    # Map sources to paths
    source_to_path = {
        "UniProtKB/SwissProt": args.out_uniprot,
        "UniProtKB/TrEMBL": args.out_uniprot,
        "Ensembl Protein": args.out_ensembl,
        "Ensembl Gene": args.out_ensembl,
        "Ensembl Transcript": args.out_ensembl,
        "RefSeq": args.out_refseq,
        "NCBI": args.out_ncbi,
        "EntrezGene": args.out_entrez,
    }

    # Open each unique path only once
    unique_paths = {p for p in source_to_path.values() if p}
    stack = ExitStack()
    path_to_handle = {p: stack.enter_context(auto_open(p, "wt")) for p in unique_paths}

    # Map each source to the corresponding handle
    outputs = {source: path_to_handle[path] for source, path in source_to_path.items() if path}
    return outputs, stack


def dump_xrefs(conf):
    with Database(conf.db) as db:
        if conf.xref_db is not None:
            xref_db = tables.open_file(conf.xref_db)
            xref_tab = xref_db.get_node("/XRef")
        else:
            xref_tab = db.db.get_node("/XRef")
        oma_id_mapper = db.id_mapper["OMA"]
        source_enum = xref_tab.get_enum("XRefSource")

        outputs, stack = open_xref_outputs(conf)
        with stack:
            for row in xref_tab.iterrows():
                source = source_enum(row["XRefSource"])
                if source not in outputs:
                    continue
                omaid = oma_id_mapper.map_entry_nr(row["EntryNr"])
                outputs[source].write(f"{omaid}\t{row['XRefId'].decode('utf-8')}\n")

        if conf.xref_db is not None:
            xref_db.close()


def dump_go(conf):
    with Database(conf.db) as db:
        oma_id_mapper = db.id_mapper["OMA"]
        with auto_open(conf.out_go, "wt") as go_fh:
            STEP = 100_000
            go_fh.write(
                dedent(
                    f"""\
            # Mapping of OMA IDs to Gene Ontology (GO) annotations for OMA release {db.get_release_name()}
            # Format: OMA ID<tab>GO term<tab>Evidence code<tab>Reference codes (comma-separated)
            """
                )
            )
            for enr in range(1, db.nr_proteins + 1, STEP):
                agg = defaultdict(set)
                for row in db.get_gene_ontology_annotations(enr, stop=enr + STEP):
                    key = (
                        oma_id_mapper.map_entry_nr(row["EntryNr"]),
                        f"GO:{row['TermNr']:07d}",
                        row["Evidence"].decode("utf-8"),
                    )
                    agg[key].add(row["Reference"].decode("utf-8"))
                for (omaid, term, evi), refs in agg.items():
                    refs_str = ",".join(refs)
                    go_fh.write(f"{omaid}\t{term}\t{evi}\t{refs_str}\n")


def dump_species_and_taxmapping(conf):
    produce_mapping_and_download_files(conf.db, conf.tax_sqlite, conf.out_species, conf.out_tax_mapping)


def _iter_protein_entries(db, chunk_size=500_000):
    """Yield :class:`~pyoma.browser.models.ProteinEntry` for every protein in the database.

    Reads the entries table in chunks of `chunk_size` rows so only a small
    number of raw numpy rows live in Python at any one time.  Each row is
    passed directly to ProteinEntry, so the entry fields are already loaded —
    only lazily accessed properties (sequence, exons, …) touch the disk later.
    """
    prot_tab = db.db.get_node("/Protein/Entries")
    xref_tab = db.db.get_node("/XRef")
    src_enum = xref_tab.get_enum("XRefSource")
    sourceCode = src_enum["SourceID"]

    n = len(prot_tab)
    for start in range(0, n, chunk_size):
        prot_chunk = prot_tab[start : start + chunk_size]
        cond_vars = {
            "enr_min": prot_chunk["EntryNr"].min(),
            "enr_max": prot_chunk["EntryNr"].max(),
            "sourceCode": sourceCode,
        }
        prot_chunk_xrefs = {
            row["EntryNr"]: row["XRefId"].decode("utf-8")
            for row in xref_tab.where(
                f"(EntryNr >= enr_min) & (EntryNr <= enr_max) & (XRefSource == sourceCode)", condvars=cond_vars
            )
        }

        for row in prot_chunk:
            yield ProteinEntry(db, row), prot_chunk_xrefs[row["EntryNr"]]


def _write_fasta_record(fh, header, seq, line_width=80):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), line_width):
        fh.write(seq[i : i + line_width])
        fh.write("\n")


def _format_main_isoform(prot):
    return "self" if prot.is_main_isoform else prot.get_main_isoform().omaid


def _format_locus(exons):
    """Format a numpy array of locus rows as the canonical location string."""
    if len(exons) == 0:
        return "n/a"
    exons = exons[exons["Start"].argsort()]
    if exons["Strand"][0] < 0:
        parts = [f"complement({e['Start']}..{e['End']})" for e in exons[::-1]]
    else:
        parts = [f"{e['Start']}..{e['End']}" for e in exons]
    return f"join({', '.join(parts)})" if len(parts) > 1 else parts[0]


def _dump_annotations(fh, db):
    """Write protein annotations genome-by-genome to avoid per-protein HDF5 queries."""
    fh.write(
        "# Dump of protein annotations (taken from original genome source, except for\n"
        '# "MainIsoform", which is the splicing form that has the most homologous matches).\n'
        "#   Note: MainIsoform=='self' indicates that this protein is OMA's main isoform,\n"
        "#         otherwise it is an OMA ID of the main isoform.\n"
        "# Format: OMA ID\tOriginal ID\tChromosome/Scaffold\tLocation\tMainIsoform\tDescription\n"
    )

    h5 = db.db
    prot_tab = h5.get_node("/Protein/Entries")
    desc_buf = h5.get_node("/Protein/DescriptionBuffer")
    xref_tab = h5.get_node("/XRef")
    src_enum = xref_tab.get_enum("XRefSource")
    source_code = src_enum["SourceID"]

    for gs in tqdm(h5.get_node("/Genome").read(), desc="genomes"):
        code = gs["UniProtSpeciesCode"].decode()
        off = int(gs["EntryOff"])
        tot = int(gs["TotEntries"])
        if tot == 0:
            continue

        # Read all entries for this genome in one slice
        entries = prot_tab[off : off + tot]

        # Source-ID xrefs for this genome (one range query, not per-protein)
        enr_min, enr_max = int(entries["EntryNr"].min()), int(entries["EntryNr"].max())
        src_xrefs = {
            int(row["EntryNr"]): row["XRefId"].decode()
            for row in xref_tab.where(
                "(EntryNr >= enr_min) & (EntryNr <= enr_max) & (XRefSource == source_code)",
                condvars={"enr_min": enr_min, "enr_max": enr_max, "source_code": source_code},
            )
        }

        # Locus table: one read for the entire genome, grouped by entry_nr
        try:
            locus_tab = h5.get_node(f"/Protein/Locus/{code}")
            locus_rows = locus_tab[:]
        except tables.NoSuchNodeError:
            locus_rows = None

        locus_by_enr = defaultdict(list)
        if locus_rows is not None:
            for row in locus_rows:
                locus_by_enr[int(row["EntryNr"])].append(row)

        for i, e in enumerate(entries):
            enr = int(e["EntryNr"])
            omaid = f"{code}{i + 1:05d}"
            src_xref = src_xrefs.get(enr, "")
            chrom = e["Chromosome"].decode()
            desc = (
                desc_buf[int(e["DescriptionOffset"]) : int(e["DescriptionOffset"]) + int(e["DescriptionLength"])]
                .tobytes()
                .decode()
            )

            exon_rows = locus_by_enr.get(enr)
            if exon_rows:
                locus_str = _format_locus(np.array(exon_rows, dtype=locus_rows.dtype))
            else:
                locus_str = "n/a"

            alt = int(e["AltSpliceVariant"])
            if alt == 0 or alt == enr:
                main_iso = "self"
            else:
                # alt is the entry_nr of the main isoform — compute its OMA ID using genome offset
                main_iso = f"{code}{alt - off:05d}"

            fh.write(f"{omaid}\t{src_xref}\t{chrom}\t{locus_str}\t{main_iso}\t{desc}\n")


def dump_sequences_and_protein_annotations(conf):
    with Database(conf.db) as db:
        with ExitStack() as stack:
            prot_fh = stack.enter_context(auto_open(conf.out_proteins, "wt")) if conf.out_proteins else None
            cdna_fh = stack.enter_context(auto_open(conf.out_cdna, "wt")) if conf.out_cdna else None
            annot_fh = stack.enter_context(auto_open(conf.out_annotations, "wt")) if conf.out_annotations else None

            if annot_fh:
                logger.info("writing protein annotations to %s", conf.out_annotations)
                _dump_annotations(annot_fh, db)
            if prot_fh:
                logger.info("writing protein sequences to %s", conf.out_proteins)
            if cdna_fh:
                logger.info("writing cDNA sequences to %s", conf.out_cdna)

            if prot_fh is not None or cdna_fh is not None:
                for prot, src_xref in tqdm(
                    _iter_protein_entries(db),
                    desc="Processing proteins",
                    total=db.db.get_node("/Protein/Entries").nrows,
                ):
                    if prot_fh is not None:
                        seq = prot.sequence
                        if seq:
                            _write_fasta_record(prot_fh, src_xref, seq)
                    if cdna_fh is not None:
                        cdna = prot.cdna
                        if cdna:
                            _write_fasta_record(cdna_fh, src_xref, cdna)


def dump_ncbi_linkout(conf):
    if not conf.skip_find:
        ncbi_linkout.prepare_linkout_files(conf.out, conf.db)
    if conf.push:
        ncbi_linkout.copy_to_ncbi(conf.out, conf.password)


def dump_uniprot_crosslinks(conf):
    with auto_open(conf.out, "wt", encoding="utf-8", newline="") as fout:
        uniprot_linkout.write_omagrp_xref_map(fout, conf.db)


def _load_oma_group_data(db):
    groups = collections.defaultdict(list)
    for row in tqdm(db.get_hdf5_handle().get_node("/Protein/Entries").where("OmaGroup > 0"), desc="Loading OMA groups"):
        groups[int(row["OmaGroup"])].append(int(row["EntryNr"]))
    meta_tab = db.get_hdf5_handle().get_node("/OmaGroups/MetaData")
    kw_buf = db.get_hdf5_handle().get_node("/OmaGroups/KeywordBuffer")
    fingerprints = {}
    keywords = {}
    for row in tqdm(meta_tab, desc="Loading group metadata"):
        grp_nr = int(row["GroupNr"])
        fingerprints[grp_nr] = row["Fingerprint"].decode()
        kw = (
            kw_buf[int(row["KeywordOffset"]) : int(row["KeywordOffset"]) + int(row["KeywordLength"])].tobytes().decode()
        )
        if kw:
            keywords[grp_nr] = kw
    xref_tab = db.get_hdf5_handle().get_node("/XRef")
    source_code = xref_tab.get_enum("XRefSource")["SourceID"]
    gene_ids = {
        int(row["EntryNr"]): row["XRefId"].decode()
        for row in tqdm(
            xref_tab.where("XRefSource == source_code", condvars={"source_code": source_code}), desc="Loading gene IDs"
        )
    }
    nr_species = len(db.get_hdf5_handle().get_node("/Genome"))
    return groups, fingerprints, keywords, gene_ids, nr_species


def _dump_oma_groups_txt(out, db, groups, fingerprints, gene_ids, nr_species):
    with auto_open(out, "wt") as fh:
        fh.write(
            f"# Orthologous groups from OMA release of {db.get_release_name()}\n"
            f"# This release has {len(groups)} groups covering {sum(len(z) for z in groups.values())} proteins from {nr_species} species\n"
            "# Format: group number<tab>Fingerprint<tab>tab-separated list of OMA Entry IDs\n"
        )
        id_mapper = db.id_mapper["OMA"]
        for group_nr in range(1, len(groups) + 1):
            oma_ids = groups[group_nr]
            members = "\t".join(map(id_mapper.map_entry_nr, oma_ids))
            fh.write(f"{group_nr}\t{fingerprints.get(group_nr, 'n/a')}\t{members}\n")


def _dump_oma_groups_orthoxml(out, db, groups, fingerprints, gene_ids, nr_species):
    from html import escape as xml_escape

    nr_proteins_total = db.get_hdf5_handle().get_node("/Protein/Entries").nrows
    with auto_open(out, "wt") as fh:
        fh.write(
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3"'
            f' origin="OMA" originVersion="{xml_escape(db.get_release_name())}">\n'
            f" <notes>\n"
            f'  <stats nrGroups="{len(groups)}" nrSpecies="{nr_species}"'
            f' nrProt="{nr_proteins_total}" nrProtInGroups="{sum(len(z) for z in groups.values())}" />\n'
            f" </notes>\n"
        )

        logger.info("Writing species sections for %d genomes...", nr_species)
        for genome in tqdm(db.get_hdf5_handle().get_node("/Genome"), desc="Writing genome sections"):
            entry_off = int(genome["EntryOff"])
            nr_entries = int(genome["TotEntries"])
            sci_name = genome["SciName"].decode()
            ncbi_taxid = int(genome["NCBITaxonId"])
            db_version = genome["Release"].decode()
            sp_code = genome["UniProtSpeciesCode"].decode()
            enr_min = entry_off + 1
            enr_max = entry_off + nr_entries

            buf = [
                f' <species name="{xml_escape(sci_name)}" NCBITaxId="{ncbi_taxid}">\n',
                f'  <database name="{xml_escape(sci_name)}" version="{xml_escape(db_version)}">\n',
                "   <genes>\n",
            ]
            for enr in range(enr_min, enr_max + 1):
                oma_id = f"{sp_code}{enr - entry_off:05d}"
                gene_id = xml_escape(gene_ids.get(enr, ""))
                buf.append(f'    <gene id="{enr}" ')
                if gene_id != "":
                    gene_id = gene_id.split("{")[0].strip()
                    buf.append(f'geneId="{gene_id}" ')
                buf.append(f'protId="{oma_id}"/>\n')
            buf.extend(["   </genes>\n", "  </database>\n", " </species>\n"])
            fh.writelines(buf)

        logger.info("Writing %d OMA groups...", len(groups))
        fh.write(" <groups>\n")
        buf = []
        for grp_nr in tqdm(range(1, len(groups) + 1), desc="Writing groups"):
            buf.append(f'  <orthologGroup id="{grp_nr}">\n')
            for j in groups[grp_nr]:
                buf.append(f'   <geneRef id="{j}" />\n')
            buf.append(f'   <notes><fingerprint id="{xml_escape(fingerprints[grp_nr])}" /></notes>\n')
            buf.append("  </orthologGroup>\n")
            if len(buf) >= 50_000:
                fh.writelines(buf)
                buf.clear()
        if buf:
            fh.writelines(buf)
        fh.write(" </groups>\n")
        fh.write("</orthoXML>\n")


def _dump_oma_groups_keywords(out_keywords, nr_groups, keywords):
    with auto_open(out_keywords, "wt") as fh:
        for grp_nr in range(1, nr_groups + 1):
            fh.write(f"{grp_nr}\t{keywords.get(grp_nr, '-')}\n")


def dump_oma_groups(conf):
    with Database(conf.db) as db:
        groups, fingerprints, keywords, gene_ids, nr_species = _load_oma_group_data(db)
        if conf.out_orthoxml is not None:
            _dump_oma_groups_orthoxml(conf.out_orthoxml, db, groups, fingerprints, gene_ids, nr_species)
        if conf.out_txt is not None:
            _dump_oma_groups_txt(conf.out_txt, db, groups, fingerprints, gene_ids, nr_species)
        if conf.out_keywords is not None:
            _dump_oma_groups_keywords(conf.out_keywords, len(groups), keywords)


def _write_roothog_descriptions(fh, h5, release_char=""):
    """Write the keywords of all root HOGs as ``RootHOG<tab>Description`` lines.

    Reads /RootHOG/MetaData and /RootHOG/KeywordBuffer in one go each rather than
    calling :meth:`Database.get_roothog_keywords` per family.
    """
    meta = h5.get_node("/RootHOG/MetaData").read()
    kw_buf = h5.get_node("/RootHOG/KeywordBuffer").read().tobytes()
    meta = meta[np.argsort(meta["FamNr"], kind="stable")]
    n = 0
    for row in meta:
        off = int(row["KeywordOffset"])
        kw = kw_buf[off : off + int(row["KeywordLength"])].decode()
        fh.write(f"HOG:{release_char}{int(row['FamNr']):07d}\t{kw}\n")
        n += 1
    return n


def dump_roothog_descriptions(conf):
    with Database(conf.db) as db:
        with auto_open(conf.out, "wt") as fh:
            fh.write(
                f"# Descriptions of root HOGs from OMA release {db.get_release_name()}\n"
                "# Format: RootHOG<tab>Description\n"
            )
            n = _write_roothog_descriptions(fh, db.get_hdf5_handle(), db.release_char)
        logger.info("wrote %d root HOG descriptions to %s", n, conf.out)


def dump_id_history(conf):
    """Write IDHistory-{old}-to-{new}.txt.gz files for each previous-release database."""
    os.makedirs(conf.out, exist_ok=True)

    with tables.open_file(conf.db, "r") as new_h5:
        new_release = str(new_h5.get_node_attr("/", "oma_version"))
        new_genome_tab = new_h5.get_node("/Genome").read()
        new_prot_tab = new_h5.get_node("/Protein/Entries")

        new_genomes = {row["UniProtSpeciesCode"]: row for row in new_genome_tab}

        for old_db_path in conf.prev_dbs:
            with tables.open_file(old_db_path, "r") as old_h5:
                old_release = str(old_h5.get_node_attr("/", "oma_version"))
                old_genome_tab = old_h5.get_node("/Genome").read()
                old_prot_tab = old_h5.get_node("/Protein/Entries")

                out_fn = os.path.join(conf.out, f"IDHistory-{old_release}-to-{new_release}.txt.gz")
                logger.info("writing %s", out_fn)

                with auto_open(out_fn, "wt") as fh:
                    fh.write(
                        f"# OMA IDs mapping from release {old_release} to {new_release} of updated genomes.\n"
                        "#\n"
                        "#  The OMA IDs of an updated genome are not stable. We provide here a\n"
                        "#  mapping of those proteins, which still have an identical protein sequence. We\n"
                        "#  do not provide a mapping for proteins which changed in their sequence.\n"
                        "#  Species which were not updated are not included in this mapping.\n"
                        "# Format: old OMA ID<tab>new OMA ID\n"
                    )

                    n_mapped = 0
                    for old_gs in old_genome_tab:
                        code = old_gs["UniProtSpeciesCode"]
                        if code not in new_genomes:
                            continue
                        new_gs = new_genomes[code]
                        if old_gs["Date"] == new_gs["Date"]:
                            continue

                        new_off = int(new_gs["EntryOff"])
                        new_tot = int(new_gs["TotEntries"])
                        new_prots = new_prot_tab[new_off : new_off + new_tot]
                        new_md5_to_pnr = {}
                        for i, row in enumerate(new_prots):
                            md5 = bytes(row["MD5ProteinHash"])
                            if md5 not in new_md5_to_pnr:
                                new_md5_to_pnr[md5] = i + 1

                        old_off = int(old_gs["EntryOff"])
                        old_tot = int(old_gs["TotEntries"])
                        old_prots = old_prot_tab[old_off : old_off + old_tot]
                        code_str = code.decode()
                        for i, row in enumerate(old_prots):
                            pnr_new = new_md5_to_pnr.get(bytes(row["MD5ProteinHash"]))
                            if pnr_new is not None:
                                fh.write(f"{code_str}{i + 1:05d}\t{code_str}{pnr_new:05d}\n")
                                n_mapped += 1

                logger.info("wrote %d ID mappings to %s", n_mapped, out_fn)


def _load_pairwise_ortholog_lookup(db):
    """Build EntryNr-indexed lookup arrays (OMA ID, OmaHOG, OmaGroup) for all main isoforms.

    VPairs only ever reference main isoforms, so restricting the scan to them keeps the
    lookup as small as possible while still covering every EntryNr that can appear in a
    pairwise relation.
    """
    entries_tab = db.get_hdf5_handle().get_node("/Protein/Entries")
    id_mapper = db.id_mapper["OMA"]
    n = entries_tab.nrows

    oma_id = [""] * (n + 1)
    oma_hog = np.empty(n + 1, dtype=object)
    oma_hog[:] = b""
    oma_group = np.zeros(n + 1, dtype=np.uint32)

    for row in tqdm(
        entries_tab.where("(AltSpliceVariant == EntryNr) | (AltSpliceVariant == 0)"),
        desc="Loading main isoforms",
    ):
        enr = int(row["EntryNr"])
        oma_id[enr] = id_mapper.map_entry_nr(enr)
        oma_hog[enr] = bytes(row["OmaHOG"])
        oma_group[enr] = row["OmaGroup"]
    # OMA IDs are near-fixed-width (species code + zero-padded index), but the index
    # is not truncated for genomes with >=100'000 entries, so let numpy infer the
    # true max width here rather than assuming one - a hardcoded fixed dtype would
    # silently truncate any longer id. OMA IDs are ASCII, so bytes (S) is more
    # compact than unicode (U); consumers must decode() the resulting elements.
    oma_id = np.array(oma_id, dtype=bytes)
    return oma_id, oma_hog, oma_group


def _iter_genome_codes_with_vpairs(db):
    """Yield (code, nrows) for every genome that has a VPairs table.

    Deliberately does *not* hand back the open Table node: reading a table's shape
    is metadata-only, but any actual row read populates a per-dataset HDF5 chunk
    cache that is only released once the node is closed.
    """
    h5 = db.get_hdf5_handle()
    for gs in h5.get_node("/Genome").read():
        code = gs["UniProtSpeciesCode"].decode()
        try:
            vp_tab = h5.get_node(f"/PairwiseRelation/{code}/VPairs")
        except tables.NoSuchNodeError:
            continue
        nrows = vp_tab.nrows
        vp_tab.close()
        yield code, nrows


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def dump_pairwise_orthologs(conf):
    """Write all pairwise orthologs (VPairs) as a TSV file.

    Format: OMA_ID1\tOMA_ID2\tRelType\tOMA_Group\tOMA_HOG

    OMA_Group is only reported if both entries are members of the same (non-zero) OMA Group.
    OMA_HOG is the deepest common sub-HOG prefix, if the two entries qualify as orthologous
    according to :func:`pyoma.browser.hoghelper.are_orthologous`. Cells that do not qualify
    are reported as "n/a".
    """
    chunk_size = 2_000_000
    logger.debug("RSS before Database(): %.1f MB", rss_mb())
    with Database(conf.db) as db:
        logger.debug("RSS after Database(): %.1f MB", rss_mb())
        logger.info("Loading OMA-ID / OmaHOG / OmaGroup lookup for main isoforms...")
        oma_id, oma_hog, oma_group = _load_pairwise_ortholog_lookup(db)
        logger.debug("RSS after lookup build: %.1f MB", rss_mb())

        genome_vpair_info = list(_iter_genome_codes_with_vpairs(db))
        logger.debug("RSS after genomes_vtab_list: %.1f MB", rss_mb())
        n_proteins = db.get_hdf5_handle().get_node("/Protein/Entries").nrows
        n_pairs = sum(nrows for _, nrows in genome_vpair_info) // 2
        n_species = len(genome_vpair_info)

        reltype_map = None
        n_written = 0
        with auto_open(conf.out, "wt") as fh:
            fh.write(
                f"# All orthologous pairs from OMA release {db.get_release_name()}\n"
                f"# This release has {n_pairs} pairs covering {n_proteins} proteins from {n_species} species\n"
                "# Format: Protein 1<tab>Protein 2<tab>Orthology type<tab>OMA group (if any)<tab>OMA HOG (if any)\n"
                "#\n"
                "#   OMA_Group is only reported if both proteins are members of the same OMA Group (SOG).\n"
                "#   OMA_HOG is the deepest common sub-HOG prefix, if the two proteins belong to the same HOG.\n"
                "#   OMA_Group and OMA_HOG are reported as 'n/a' if the two proteins do not belong to "
                "#   the same group or HOG, respectively.\n"
                "# Every pair is listed only once, and in no particular order.\n"
            )
            h5 = db.get_hdf5_handle()
            with tqdm(total=n_pairs, desc="nr of pairs") as pbar:
                for i, (code, nrows) in enumerate(genome_vpair_info):
                    if i % 100 == 0:
                        logger.debug("RSS before processing %s [%d] vtab: %.1f MB", code, i, rss_mb())
                    vp_tab = h5.get_node(f"/PairwiseRelation/{code}/VPairs")
                    try:
                        if reltype_map is None:
                            type_enum = vp_tab.get_enum("RelType")
                        reltype_map = {val: name for name, val in type_enum}

                        for start in range(0, nrows, chunk_size):
                            chunk = vp_tab[start : start + chunk_size]
                            mask = chunk["EntryNr1"] > chunk["EntryNr2"]
                            sub = chunk[mask]
                            if len(sub) == 0:
                                continue

                            e1 = sub["EntryNr1"]
                            e2 = sub["EntryNr2"]
                            ids1 = oma_id[e1]
                            ids2 = oma_id[e2]
                            hogs1 = oma_hog[e1]
                            hogs2 = oma_hog[e2]
                            g1 = oma_group[e1]
                            g2 = oma_group[e2]
                            group_match = (g1 == g2) & (g1 > 0)
                            rel_types = sub["RelType"]

                            buf = []
                            for i in range(len(sub)):
                                ortho_type = reltype_map[int(rel_types[i])]
                                group_str = str(int(g1[i])) if group_match[i] else "n/a"

                                h1, h2 = hogs1[i], hogs2[i]
                                hog_str = "n/a"
                                if h1 and h2:
                                    if h1 == h2:
                                        hog_str = h1.decode()
                                    else:
                                        common = are_orthologous(h1, h2, return_group=True)
                                        if common is not False:
                                            hog_str = common

                                buf.append(
                                    f"{ids1[i].decode()}\t{ids2[i].decode()}\t{ortho_type}\t{group_str}\t{hog_str}\n"
                                )
                            fh.writelines(buf)
                            n_written += len(sub)
                            pbar.update(len(sub))
                    finally:
                        vp_tab.close()
        logger.info("wrote %d pairwise ortholog relations to %s", n_written, conf.out)
        logger.info("RSS after writing pairwise orthologs: %.1f MB", rss_mb())


def parse_command_line_args():
    parser = ArgumentParser(description="Dump various files from the OMA database")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity")

    subparsers = parser.add_subparsers(title="Commands")

    xrefs_parser = subparsers.add_parser(
        "xrefs", help="Dump cross-reference mappings", formatter_class=ArgumentDefaultsHelpFormatter
    )
    xrefs_parser.set_defaults(func=dump_xrefs)
    xrefs_parser.add_argument("--db", required=True, help="Path to database")
    xrefs_parser.add_argument(
        "--xref-db", required=False, help="Path to xref database. If not provided, the main database is used."
    )
    xrefs_parser.add_argument("--out-uniprot", help="Output path for the uniprot cross-reference mappings (if enabled)")
    xrefs_parser.add_argument("--out-ensembl", help="Output path for the ensembl cross-reference mappings (if enabled)")
    xrefs_parser.add_argument("--out-refseq", help="Output path for the refseq cross-reference mappings (if enabled)")
    xrefs_parser.add_argument("--out-ncbi", help="Output path for the ncbi cross-reference mappings (if enabled)")
    xrefs_parser.add_argument("--out-entrez", help="Output path for the entrez cross-reference mappings (if enabled)")

    go_parser = subparsers.add_parser("go", help="Dump GO annotations", formatter_class=ArgumentDefaultsHelpFormatter)
    go_parser.set_defaults(func=dump_go)
    go_parser.add_argument("--db", required=True, help="Path to database")
    go_parser.add_argument("--out-go", help="Output path for the GO annotations")

    species_parser = subparsers.add_parser(
        "species", help="Dump species information", formatter_class=ArgumentDefaultsHelpFormatter
    )
    species_parser.set_defaults(func=dump_species_and_taxmapping)
    species_parser.add_argument("--db", required=True, help="Path to database")
    species_parser.add_argument("--tax-sqlite", required=True, help="Path to the taxonomy.sqlite file")
    species_parser.add_argument(
        "--out-species", required=True, help="Output path for the species information in tsv format"
    )
    species_parser.add_argument(
        "--out-tax-mapping", required=True, help="Path where the GTDB / NCBI taxonomy mapping is stored as pickle file."
    )

    seq_parser = subparsers.add_parser(
        "sequences",
        help="Dump protein/cDNA sequences as FASTA and protein annotations as TSV",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    seq_parser.set_defaults(func=dump_sequences_and_protein_annotations)
    seq_parser.add_argument("--db", required=True, help="Path to database")
    seq_parser.add_argument("--out-proteins", help="Output path for protein sequences in FASTA format")
    seq_parser.add_argument("--out-cdna", help="Output path for cDNA sequences in FASTA format")
    seq_parser.add_argument("--out-annotations", help="Output path for protein annotations in TSV format")

    ncbi_parser = subparsers.add_parser(
        "ncbi-linkout",
        help="Generate and optionally push NCBI LinkOut XML files",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    ncbi_parser.set_defaults(func=dump_ncbi_linkout)
    ncbi_parser.add_argument("--db", required=True, help="Path to the HDF5 database")
    ncbi_parser.add_argument("--out", default="./", help="Folder where output XML files are written")
    ncbi_parser.add_argument("--skip-find", action="store_true", help="Skip generating the linkout files")
    ncbi_parser.add_argument("--push", action="store_true", help="Push XML files to NCBI FTP host")
    ncbi_parser.add_argument("-p", "--password", help="Password for the NCBI FTP account")

    uniprot_parser = subparsers.add_parser(
        "uniprot-crosslinks",
        help="Generate UniProtKB → OMA group crosslink mapping",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    uniprot_parser.set_defaults(func=dump_uniprot_crosslinks)
    uniprot_parser.add_argument("--db", required=True, help="Path to the HDF5 database")
    uniprot_parser.add_argument("--out", default="./UniProtKB-OMA.txt.gz", help="Output file for the crosslink mapping")

    groups_parser = subparsers.add_parser(
        "oma-groups", help="Dump OMA groups in TXT or OrthoXML format", formatter_class=ArgumentDefaultsHelpFormatter
    )
    groups_parser.set_defaults(func=dump_oma_groups)
    groups_parser.add_argument("--db", required=True, help="Path to database")
    groups_parser.add_argument("--out-orthoxml", help="Output path for OMA groups in OrthoXML format")
    groups_parser.add_argument("--out-txt", help="Output path for OMA groups in txt format")
    groups_parser.add_argument("--out-keywords", help="Output path for OMA group keywords (TSV: group_nr<tab>keywords)")

    roothog_desc_parser = subparsers.add_parser(
        "roothog-descriptions",
        help="Dump the descriptions (keywords) of root HOGs as TSV (RootHOG<tab>Description)",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    roothog_desc_parser.set_defaults(func=dump_roothog_descriptions)
    roothog_desc_parser.add_argument("--db", required=True, help="Path to database")
    roothog_desc_parser.add_argument(
        "--out", default="roothog-descriptions.txt.gz", help="Output path for the root HOG descriptions"
    )

    pairwise_orthologs_parser = subparsers.add_parser(
        "pairwise-orthologs",
        help="Dump pairwise orthologs (VPairs) as a TSV file, together with OMA Group and OmaHOG information",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    pairwise_orthologs_parser.set_defaults(func=dump_pairwise_orthologs)
    pairwise_orthologs_parser.add_argument("--db", required=True, help="Path to database")
    pairwise_orthologs_parser.add_argument(
        "--out", required=True, help="Output path for the pairwise orthologs TSV file"
    )

    id_history_parser = subparsers.add_parser(
        "id-history",
        help="Produce IDHistory files mapping OMA IDs across updated genomes between releases",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    id_history_parser.set_defaults(func=dump_id_history)
    id_history_parser.add_argument("--db", required=True, help="Path to the current (new) HDF5 database")
    id_history_parser.add_argument(
        "--prev-dbs",
        nargs="+",
        required=True,
        metavar="PREV_DB",
        help="Paths to previous-release HDF5 databases",
    )
    id_history_parser.add_argument(
        "--out", default=".", help="Output directory for IDHistory-{old}-to-{new}.txt.gz files"
    )

    conf = parser.parse_args()
    if not hasattr(conf, "func"):
        parser.print_usage()
        sys.exit(1)
    return conf


def dump_files():
    conf = parse_command_line_args()
    setup_logging(conf)
    logger.info("Command line options: %s", str(conf))
    if not sys.warnoptions and not getattr(conf, "verbose", 0) >= 1:
        warnings.simplefilter("ignore", category=PerformanceWarning)
        warnings.simplefilter("ignore", category=RuntimeWarning)
    conf.func(conf)


if __name__ == "__main__":
    dump_files()
