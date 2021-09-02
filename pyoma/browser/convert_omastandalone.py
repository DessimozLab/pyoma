from .convert import *
from .compute_cache import compute_and_store_cached_data
import os


class StandaloneExporter(DarwinExporter):
    DRW_CONVERT_FILE = os.path.abspath(os.path.splitext(__file__)[0] + ".drw")

    def __init__(self, root, name, **kwargs):
        os.environ["DARWIN_BROWSERDATA_PATH"] = os.path.abspath(root)
        super(StandaloneExporter, self).__init__(name, **kwargs)
        self.transformed = False
        self.cache_dir = os.path.join(os.getenv("DARWIN_BROWSERDATA_PATH"), "pyoma")

    def add_homologs(self):
        self.assert_cached_results()
        for gs in self.h5.root.Genome.iterrows():
            genome = gs["UniProtSpeciesCode"].decode()
            rel_node_for_genome = self._get_or_create_node(
                "/PairwiseRelation/{}".format(genome)
            )
            if "homologs" not in rel_node_for_genome:
                pass

    def get_version(self):
        # TODO: obtain real version
        return "OmaStandalone; 1.0.x"

    def assert_cached_results(self):
        if not self.transformed:
            res = self.call_darwin_export(
                "TransformDataToCache('{}');".format(self.cache_dir)
            )
            if res != "success":
                raise DarwinException("could not transform data from darwin", "")
            self.transformed = True
            os.environ["DARWIN_NETWORK_SCRATCH_PATH"] = os.getenv(
                "DARWIN_BROWSERDATA_PATH"
            )
            common.package_logger.info("successfully transformed data to json")

    def add_orthologs(self):
        self.assert_cached_results()
        for gs in self.h5.root.Genome.iterrows():
            genome = gs["UniProtSpeciesCode"].decode()
            rel_node_for_genome = self._get_or_create_node(
                "/PairwiseRelation/{}".format(genome)
            )
            if "VPairs" not in rel_node_for_genome:
                cache_file = os.path.join(
                    os.getenv("DARWIN_NETWORK_SCRATCH_PATH", ""),
                    "pyoma",
                    "vps",
                    "{}.txt.gz".format(genome),
                )
                if os.path.exists(cache_file):
                    data = load_tsv_to_numpy((cache_file, 0, 0, False))
                else:
                    # fallback to read from VPsDB
                    data = self.call_darwin_export("GetVPsForGenome({})".format(genome))

                vp_tab = self.h5.create_table(
                    rel_node_for_genome,
                    "VPairs",
                    tablefmt.PairwiseRelationTable,
                    expectedrows=len(data),
                )
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, vp_tab)
                self._write_to_table(vp_tab, data)
                vp_tab.cols.EntryNr1.create_csindex()

    def add_hogs(self, **kwargs):
        fn = "HierarchicalGroups.orthoxml"
        hog_file = os.path.join(os.environ["DARWIN_BROWSERDATA_PATH"], "Output", fn)
        hog_cache_dir = os.path.join(self.cache_dir, "split_hogs")
        for tree_file in (
            "ManualSpeciesTree.nwk",
            "EstimatedSpeciesTree.nwk",
            "LineageSpeciesTree.nwk",
        ):
            tree_filename = os.path.join(
                os.environ["DARWIN_BROWSERDATA_PATH"], "Output", tree_file
            )
            if os.path.exists(tree_filename):
                self.logger.info("Use " + tree_filename + " as HOG backbone tree file")
                break
        hog_treefile = None
        if os.path.exists(tree_filename):
            hog_treefile = tree_filename
        return super().add_hogs(
            hog_path=hog_cache_dir, hog_file=hog_file, tree_filename=hog_treefile
        )

    def _get_genome_database_paths(self):
        return self.call_darwin_export("GetGenomeFileNames();")

    def xref_databases(self):
        return self._get_genome_database_paths()


def mark_isoforms(dbfn):
    from .db import Database, OmaIdMapper

    used_splice_fn = os.path.join(
        os.getenv("DARWIN_BROWSERDATA_PATH"), "Output", "used_splicing_variants.txt"
    )
    if not os.path.exists(used_splice_fn):
        common.package_logger.info(
            "no splicing output file found. Assume not splice variants"
        )
        return

    db = Database(dbfn)
    idmapper = OmaIdMapper(db)
    main_variants = {}
    with open(used_splice_fn) as fh:
        for line in fh:
            try:
                main_variant = line.split("\t")[1].split("|")[0].strip()
                common.package_logger.debug(main_variant)
                main_variants[main_variant] = idmapper.omaid_to_entry_nr(main_variant)
            except Exception:
                common.package_logger.warning("cannot convert line: {}".format(line))
                pass
    common.package_logger.info(
        "found {} main splicing variants".format(len(main_variants))
    )
    common.package_logger.debug(main_variants)
    splice_column = (
        db.get_hdf5_handle().get_node("/Protein/Entries").col("AltSpliceVariant")
    )
    for file in os.scandir(os.path.join(os.getenv("DARWIN_BROWSERDATA_PATH"), "DB")):
        if not file.name.endswith(".splice"):
            continue
        common.package_logger.warning("handling {}".format(file.name))
        with open(file) as fh:
            for line in fh:
                splice_variants = [z.strip() for z in line.split(";")]
                main = [v for v in splice_variants if v in main_variants]
                if len(main) != 1:
                    common.package_logger.warning(
                        "not a single main variant for {}: {}".format(
                            splice_variants, main
                        )
                    )
                    continue
                for v in splice_variants:
                    enr = idmapper.omaid_to_entry_nr(v)
                    splice_column[enr] = main_variants[main[0]]
    db.close()
    with tables.open_file(dbfn, "a") as h5:
        tab = h5.get_node("/Protein/Entries")
        tab.modify_column(column=splice_column, colname="AltSpliceVariant")


def import_oma_run(path, outfile, domains=None, log_level="INFO"):
    log = getLogger(log_level)
    x = StandaloneExporter(path, outfile, logger=log, mode="write")
    x.add_version()
    x.add_species_data()
    x.add_orthologs()
    x.add_proteins()
    x.add_hogs()
    x.add_xrefs()
    if domains is None:
        domains = ["file://dev/null"]
    else:
        domains = list(
            map(lambda url: "file://" + url if url.startswith("/") else url, domains)
        )
    log.info("loading domain annotations from {}".format(domains))
    x.add_domain_info(
        filter_duplicated_domains(
            only_pfam_or_cath_domains(
                itertools.chain.from_iterable(map(iter_domains, domains))
            )
        )
    )
    x.add_domainname_info(
        itertools.chain(
            CathDomainNameParser(
                "http://download.cathdb.info/cath/releases/latest-release/"
                "cath-classification-data/cath-names.txt"
            ).parse(),
            PfamDomainNameParser(
                "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
            ).parse(),
        )
    )
    x.add_canonical_id()
    x.add_group_metadata()
    x.add_hog_domain_prevalence()
    x.add_roothog_metadata()
    x.close()

    x = StandaloneExporter(path, outfile, logger=log)
    x.create_indexes()
    x.add_sequence_suffix_array()
    x.update_summary_stats()
    x.add_per_species_aux_groupdata()
    x.close()
    mark_isoforms(os.path.join(path, outfile))

    compute_and_store_cached_data(x.h5.filename, nr_procs=min(os.cpu_count(), 12))


if __name__ == "__main__":
    import_oma_run("~/Repositories/OmaStandalone", "oma.h5")
