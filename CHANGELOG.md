# Changelog

All notable changes to pyoma are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and pyoma uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.15.0] - 2026-10-30

### Added

- New `oma-dump` subcommands: `roothog-descriptions` (root HOG keywords as
  TSV), `pairwise-orthologs` (VPairs as TSV, with OMA Group and OMA HOG),
  `id-history` (IDHistory files mapping OMA IDs of updated genomes between
  releases), `ncbi-linkout` and `uniprot-crosslinks`.
- `oma-dump oma-groups --out-keywords` writes the OMA group keywords as TSV.
- `oma-dump sequences` writes OMA IDs in the protein and cDNA FASTA headers.
- `oma-dump xrefs --out-source` writes a mapping from OMA IDs to the source
  IDs and accessions of each protein.
- `Database.iter_oma_groups()` to iterate over all OMA groups.
- `StructureDB.iter_seq_3di()` to read 3Di sequences efficiently in bulk.
- Large new test suites for the database build (`builder.py`,
  `hogconvert.py`), the dumps and the taxonomy tree consumers.

### Changed

- Requires numpy >= 2.0. Supported Python versions are now 3.10 to 3.14
  (3.9 is dropped, 3.13 and 3.14 are new).
- ete3 is replaced by ete4. The `create_db` extra now requires `ete4>=4.3`,
  `pyham>=2.0` and `omataxonomy>=0.5`. HOG profiles can differ slightly for
  HOGs where pyham 2 counts genes, duplications and losses differently.
- The species tree of the HOG profiles is stored as a table instead of a
  pickled ete3 tree, so importing `pyoma.browser.db` no longer needs ete3.
  Files in the old format are converted on load. Profiles built with this
  version cannot be read by older pyoma versions.
- The HOG ID mapping between releases (`hogidmap`) is reworked: it uses an
  LSH forest, can compute exact Jaccard similarities, shuffles families of
  unknown size and stops workers that exceed a timeout.
- `CanonicalId` uses the SwissProt entry name where available.
- NCBI and UniProt linkout generation moved into the `pyoma.browser.linkout`
  package.
- The version is derived from git tags with `setuptools_scm`.
- All dependencies have lower bounds, and CI also tests the minimal
  versions (on Python 3.10) besides the locked and newest versions. CI uses
  uv; `tox.ini` is removed.
- Faster collection of GO summary and per-species metadata, and faster
  dumping of protein annotations.

### Fixed

- Exon order for genes on the negative strand in the protein annotation
  dump.
- numpy uint64/float64 upcast when slicing in `get_sequence` and `get_cdna`.
- Name of the `/XRefIndex/XRef_reduced` table.
- Taxonomy mapping when genomes with negative GenomeIds exist but no GTDB
  data.
- UniProt-OMA linkout only uses accessions from SwissProt and TrEMBL.
- `PerGenomeOMAGroupAuxData` for databases without OMA groups, and importing
  OMA groups from FastOMA runs.
- Jaccard computation with duplicated protein IDs in `hogidmap`.
- Memory use of the pairwise orthologs dump.
- Database build: `OmaGroupsProvider` without OMA groups, close paralogs of
  genome pairs without shared paralogs, choice of the main splice variant
  from pairwise orthologs, root HOGs in `HOGtoHDF5.process_augmented_hog`,
  and masked exceptions in `HOGtoHDF5.__exit__`.

## [0.14.2] - 2026-07-27

### Fixed

- HOG ID mapping (`hogidmap`) works with `datasketch` >= 2.0. The MinHash `scheme` attribute is set to `legacy` when it was not stored before; otherwise the scheme the MinHash was built with is kept.
- Type annotations in `idmapper` work again on Python 3.9.
- The GO enrichment test is now deterministic.

## [0.14.1] - 2026-07-27

### Added

- Initial support for protein structure data: a new `pyoma.browser.structure` module with `StructureDB` and `StructureInfo` (3Di sequence and its source). `Database.load_structure_db()` loads this data from an external HDF5 file, by default `OmaServer.structure.h5`.
- Prefix and substring search in the gene name and protein name lookups (`GeneNamesLookup` / `ProteinNamesLookup`).
- `oma-dump oma-groups` exports OMA Groups in text and OrthoXML formats.
- Synteny graphs include the orientation and an orientation score for each edge.

### Changed

- `Database.get_conversion_date()` accepts an optional `kind` argument.

### Fixed

- Handling of the `conversion_start`/`conversion_end` PyTables attributes.

## [0.14.0] - 2026-05-20

### Added

- New database build pipeline in `pyoma.browser.build`, with the `oma-build` command. It covers genome and protein import, sequence indexes, HOG conversion (including FastOMA output), domain annotations, xref fetching, filtering and mapping, Gene Ontology import, OMA Group based GO prediction, fingerprints, keywords, homoeologs and close paralogs, and the cache build. Many steps run in parallel with `mpire`.
- New `oma-dump` command that writes the main download files, such as sequences, protein annotations and genome information.
- Suffix-array based sequence search, with the sequence buffer stored in a separate memory-mapped file, plus faster k-mer approximate search. Xrefs can be mapped inexactly from k-mer counts alone, without an alignment.
- Species divergence times, a `FamIndex` structure for fast loading of family members, summary metadata, and synteny age in the `ExtantSyntenyRels` table.
- A `hogid-map` CLI.
- OrthoXML output can be gzipped and pretty-printed.

### Changed

- Database schema bumped to 3.8:
  - `Taxonomy` has a new `IsGenome` column.
  - `CanonicalId` was widened to 255 characters.
  - `XRef` has a new `Identity` column and new sources (`Alternative Protein Name`, `KEGG`, `AGR`).
  - Synteny tables have `Orientation` and `OrientationScore` columns.
- `CanonicalId` is chosen by a fixed `XRefSource` priority order.
- Species codes may be any 5 characters from `[A-Z0-9]`, so codes from external providers work.
- Protein IDs only need to be unique within a genome.
- Packaging moved from `setup.py` to `pyproject.toml`. Python 3.9–3.12 is required, with `numpy>=1.26,<2.0` and `tables>=3.9.0`.
- HDF5 files use blosc2 compression.
- The test database is downloaded from Zenodo.

### Fixed

- `taxmap` handles special genomes (extant species on the taxonomy path of another genome).
- Large HOG ID map data that was too long for a `VLArray` is stored by manual pickling.
- Several pipeline deadlocks and logging problems.

### Removed

- Legacy scripts `importdata.py`, `update_cached.py`, `compute_ancestral_synteny.py`, `integrate_edgehog.py`, `integrate_hogprop.py` and `taxonomymap.py`.
- The modules `compute_cache`, `ancestral_synteny` and `convert_omastandalone`.
- The `familyanalyzer` dependency.

## [0.13.5] - 2026-02-27

### Added

- `ExtantSyntenyRels` table definition.
- `ProteinEntry` model methods `is_orthologous_to` and `is_paralogous_to`, based on the OmaHOG ID.

### Fixed

- Xref search handles IDs that end in `.NN` (e.g. PomBase `SPBC11B10.09`). It first searches the full ID and only truncates the suffix if nothing is found.
- The Taxonomy root node is built reliably across schema versions.

## [0.13.4] - 2024-09-25

### Added

- `Taxonomy.newick()` has new `leaf`, `internal` and `quoted` arguments. Leaves can be labelled by scientific name or UniProt species code, internal nodes by name or taxid, and labels can be quoted instead of escaped.
- The extant synteny graph node data includes the gene end position.

### Changed

- PyPI releases use trusted publishing.

## [0.13.3] - 2024-07-03

### Changed

- Database schema bumped to 3.7: the `TaxonomyTable` has a new `Age` field and the ancestral synteny relations have a new `LCA_taxid` column.
- `numpy` is pinned to `<2` for PyTables compatibility.

## [0.13.2] - 2024-04-09

### Added

- `XRefIdMapper` can be created without search functionality, which is faster to instantiate.
- `PairwiseRelation` models accept `ProteinEntry` instances, for better performance.
- `map_to_closest_seq` logs timing statistics.

### Changed

- Python 3.12 is excluded from CI because `scikit-fuzzy` does not support it yet.

## [0.13.1] - 2023-11-10

### Fixed

- Outdated HOG IDs are resolved in search.
- Logging uses lazy evaluation.

## [0.13.0] - 2023-10-16

### Added

- GO enrichment analysis in the new `pyoma.application.enrichment` module, with a new `enrichment` extra. Ancestral GO annotations, GO term counts stored in the HDF5 file, and an information content property on the GO model.
- Extended search module: domain, GO and EC searchers, HOG ID parsing, and xref search limited to a taxonomic range (new `XRef_EntryNr_offset` array). Fuzzy genome search covers extant genomes and returns a similarity score. New `count_*` methods and a result size limit.
- Approximate sequence search can run without alignment, or with local or global alignment.
- Ancestral and extant synteny graph methods.
- New `pyoma.browser.idmapper` module, including `GeneNameAndMainDbIdMapper`.
- New `pyoma.browser.exceptions` module.
- Lineages can include a LUCA level, and LUCA is added as the root of the taxonomy data when needed.
- `TotGenes` column in the genome table. Helper scripts to integrate edgehog and hogprop results and to map GTDB to NCBI taxonomy.

### Changed

- Per-level HOGs moved from `/Hogs_per_Level/taxXXX` to `/AncestralGenomes/taxXXX/Hogs`. The `HogLevel` table has a new per-level index column.
- Python 3.8–3.11 is supported; Python 3.6 and 3.7 support was dropped. Requires `numpy>=1.20`, `tables>=3.8.0`, `numexpr>=2.8` and `property_manager`.
- Non-NCBI (negative) taxon IDs start at -2, so traversals no longer recurse forever.
- Performance improvements for xref search, OMA group sizes and table queries.

### Fixed

- Works with numexpr 2.8.5: query variables are passed separately.
- Removed uses of `numpy.bool`, which numpy 1.24 removed.
- `identify_genome` works with negative taxids.
- Unhandled `TooUnspecificQuery` errors in xref suffix search.

## [0.12.1] - 2022-04-08

### Added

- Function to extract HOG-induced orthologs between a pair of genomes.
- `oma_id` alias property for `omaid`.

### Fixed

- The fletcher32 filter is no longer used together with blosc, which HDF5 does not support.

## [0.12.0] - 2022-01-20

### Added

- New `pyoma.browser.search` module with an `XRefSearch` object. `SuffixIndex` has a new `count` method and a `limit` argument for `find`.
- `AncestralGenome` model. Genome models can be created from a species code or a taxid.
- `Database` method that loads all orthologs between two genomes efficiently.
- `Database` and `HogIDForwardMapper` can be used as context managers.
- `XRefIdMapper` can map a range of entry numbers. A mixin selects only the best ID per xref source.
- Sequence search can be limited to a set of entry numbers used as a taxon filter.
- `map_closest_seq` supports an OMAmer-like output format.
- The import can run as separate phases, with new `add_xrefs` and `project_functions` scripts.
- Compressed output files can be written.
- GitHub Actions test runner, including Python 3.10.

### Changed

- The sanity-check notebook moved to `notebooks/`, with a new `notebooks` extra.
- The fletcher32 filter was removed for now because of a PyTables issue.

### Fixed

- The OrthoXML index did not contain all families because a flush was missing.
- The suffix array can be built from empty data.
- Suffix search min/max bounds were wrong.
- Extra SwissProt xrefs can be added late in the import.

## [0.11.1] - 2021-09-02

### Changed

- Depends on `fuzzyset2` instead of the unmaintained `fuzzyset`.

### Fixed

- The `compute_cache` call in `convert_omastandalone.py`.

## [0.11.0] - 2021-08-27

### Added

- `approx_chromosome_length` method on the `Genome` model.
- Species-specific GO mapping.
- `genomes.json` stores the number of evolutionary events on terminal branches.
- Docstrings for the `models` module.
- GitHub workflow that releases the package to PyPI.

### Changed

- Python >= 3.6 is required.
- The `biopython<=1.77` restriction was removed.
- `oma_group_members` is faster.

### Fixed

- Root HOG parsing with a release character, including in the domain prevalence computation.
- Species disambiguation in per-level HOGs.
- Error handling in `callDarwinExport` when Darwin fails.

## [0.10.1] - 2021-05-28

### Added

- `nr_genes` property on the `Genome` model.
- `auto_open` helper for opening database files.

### Changed

- Main isoforms load faster.

### Fixed

- Instance check for a PyTables file handle in the `Database` constructor.

## [0.10.0] - 2021-05-10

### Added

- HOG IDs can contain a release character, which is stored as the `/oma_release_char` attribute in the HDF5 file. Use `format_hogid` to build HOG IDs; `OutdatedHogId` is raised when the release does not match.
- HOG ID mapping between releases based on MinHash LSH, computed in parallel.
- Search for outdated HOG IDs.
- Consistency check that pairwise orthologs have good scores and distances.

### Changed

- `pebble` is a dependency for building the HDF5 database.
- Tests are no longer included in the package.
- Non-informative lineage levels are removed from the possible duplication levels.

### Fixed

- Suffix search problems.
- Truncated IDs in `CanonicalId`; full IDs are used.
- Proteins at the end with many missing description values.
- The synteny graph reports a clear error when the query node is missing.

## [0.9.4] - 2021-01-25

### Added

- `Database.get_all_hogs_at_level()`, with an option to compare the HOGs against another level, and a `compare_levels` helper to compare HOGs between two taxonomic levels.
- The convert pipeline now caches `HogLevel` tables per level (computed in parallel), and the gain, loss, retained and duplication counts per genome are added to `genomes.json`.
- Start of HOG ID mapping between releases (`hogidmap` module, MinHash-based).
- `bin/update_xrefs` script to update only the cross-references in an existing HDF5 database.
- Scripts to dump HOG IDs within a clade (`dump-hogids-at-level.py`) and protein LOFT IDs (`dump-protein-loft-ids.py`).

### Changed

- `get_hog_induced_pairwise_paralogs()` is faster and supports slicing with `start`/`stop`.

### Fixed

- Wrong duplication and loss event counts.
- Parent HOGs that sit on a different subhog are now handled correctly.

## [0.9.3] - 2020-11-17

### Changed

- Faster loading of HOG members, and a simpler way to resolve the right HOG for a given level.

### Fixed

- The wrong subhog was returned for HOGs with many paralog clusters at the same level.

## [0.9.2] - 2020-11-05

### Added

- Helper function that returns the `seq_match` state of a cross-reference as a string.

## [0.9.1] - 2020-11-02

### Added

- `Database.get_subhogs_at_level()`.

## [0.9.0] - 2020-10-28

### Added

- HOG profiles (`hogprofile` package): computing profiles and searching for families with similar profiles through `Database.get_families_with_similar_hog_profile()`.
- Ancestral synteny graphs are stored in the HDF5 database; new `get_syntentic_hogs()` and `get_neighbor_hogs()` methods.
- New HOG navigation API: `get_hog()`, `get_parent_hogs()`, and `get_subhogs()` (now a generator of `HOG` instances, with `level` and internal-levels-only options). The `HOG` and `OmaGroup` models are extended.
- Augmented OrthoXML (with all taxonomic levels and taxids), root HOG descriptions and keywords, and cached family JSON and ortholog counts are stored in the database (`get_orthoxml()`, `get_roothog_keywords()`, `get_cached_family_json()`, `nr_ortholog_relations()`).
- Cross-references now report `seq_match` (exact, unchecked or modified). `resolve()` can flag IDs that belong to modified sequences. More cross-references are imported from UniProt/SwissProt (projected xrefs).
- New search features: `search_protein()`, search in descriptions, and faster xref search.
- Gene similarity within HOGs (`get_gene_similarities_hog()`), induced orthologs now include the relation type, and the protein model has `nr_exons`.
- Isoform information from OMA standalone runs is imported into HDF5.
- `map_to_closest_seq` script and a Python-based UniProt xref dumper.
- Start of an OrthoXML comparison module (`pyoma.orthoxml.compare`).

### Changed

- Python 2 support was dropped (`python_requires >= 3.5`); Python 3.8 is now supported.
- New dependencies: `datasketch`, `ete3` and `networkx`. `pyham` was added to the `create_db` extra, and `familyanalyzer>=0.7.3` is now required.
- Compatible with Biopython >= 1.78 and pandas > 0.25.
- The package is licensed under MPL 2.0 and now has a README and full package metadata. The code is formatted with black (a pre-commit config was added).

### Fixed

- `get_hog()` returned the wrong HOG when a level was given.
- Approximate sequence search failed with repetitive k-mers and when the match was at the end of the index.
- The duplication level is now placed on the branch leading to the level.
- Handling of LUCA in parent HOG and subhog lookups.

## [0.8.0] - 2020-03-06

### Added

- `HOG` and `OmaGroup` models, and `Database` methods to get information on a (sub)HOG, including `get_subhogs()`.
- `count_hog_members()`, and paging (`start`/`stop`/`step`) for HOG member retrieval.
- `get_hog_induced_pairwise_paralogs()`.
- Per-species metadata (most and least similar species, number of genes in groups), computed during conversion and stored in the HDF5 database.
- The `HogLevel` table has more information. The DB schema version was bumped.

### Fixed

- Python 3.7+ compatibility problem in conversion.

## [0.7.2] - 2019-11-25

### Added

- `Database.get_hog_induced_pairwise_orthologs()`.

### Changed

- Faster `hog_levels_of_fam`.
- 'X' is replaced by 'N' in cDNA/CDS sequences.
- Python 2.7 installs are pinned to `tables < 3.6.0`.

### Fixed

- Error when an exact sequence match has no GO annotations.

## [0.7.1] - 2019-08-30

### Added

- `Database.close()` method.
- The package version is available in `pyoma.__version__`. The pyoma version used to build a database is stored in it and logged when the database is opened.

### Changed

- Better thread safety with PyTables. Requires `tables >= 3.5.1`; the numpy pin was lifted.
- Supported Python versions updated (3.6 and 3.7 added, 3.3 removed).

## [0.7.0] - 2019-07-15

### Added

- Generic suffix-array indexer (`suffixsearch`), used for xref, GO term and EC annotation search.
- Approximate (fuzzy) search for genomes (`approx_search_genomes()`), and approximate xrefs stored in HDF5.
- Lookups for entries by EC number, domain ID and GO term (`entrynrs_with_ec_annotation()`, `entrynrs_with_domain_id()`, `entrynrs_with_go_annotation()`), plus `get_ec_annotations()`.
- Extraction of HOGs at specific levels or below.
- `dump_homoeologs` script.

### Changed

- HOG member functions always return numpy arrays with dtype information.
- `search_approx` returns an empty list instead of `None` when nothing is found.
- `PySAIS` is only imported when building indexes. Supports `DBRELEASE` dates like `04-May-2019`.
- Xrefs of modified sequences are no longer used as representative protein xrefs.

### Fixed

- Locus grammar on Python 2 and on systems other than macOS.
- Proteins without exon information.
- Several fixes for OMA standalone conversion.
- numpy 1.16.0 is excluded because of a "cannot set WRITABLE flag" error.

## [0.6.4] - 2018-12-19

### Added

- `OmaIdMapper` can identify a genome by its scientific name.

## [0.6.3] - 2018-12-14

### Added

- Script to map sequences to HOGs via the closest sequence (`map2hog`), with parallel execution.
- `xrefs` lazy property on the protein model.
- Taxonomy methods to get a sub-taxonomy and extend genomes.
- Genevestigator exporter and a script to dump GO annotations.

### Fixed

- Creating a `ProteinEntry` from a `numpy.integer`.
- Loading group metadata from the Darwin format.

## [0.6.2] - 2018-09-19

### Fixed

- The taxonomy API now correctly reports species that are also internal levels.
- Updated NCBI linkout URL for proteins.

## [0.6.1] - 2018-08-12

### Added

- PhyloXML serializer for the taxonomy, with optional genome information.

### Fixed

- Sequence-based OMA group identification is skipped for sequences with non-amino-acid characters, and edge cases are handled.

## [0.6.0] - 2018-08-08

### Added

- Homoeolog support: confidence scores computed during conversion, plus `get_homoeologs()`, `count_homoeologs()` and a `subgenome` attribute on proteins.
- Exon structure: `get_exons()`, `ExonStructure`/`Exon` models, and a lark-based locus parser.
- Isoform information: `get_splicing_variants()`, `is_main_isoform` and `alternative_isoforms` on the protein model.
- `all_proteins_of_genome()`, `get_gene_ontology_annotations()` for a range of entries (as DataFrame or GAF), `load_gene_ontology()` with relation subsets, and `get_induced_taxonomy()`.
- New xref types: WormBase, EnsemblGenomes, EnsemblGenes and PDB. Xrefs are deduplicated and kept in a canonical order.
- `NrMembers` in the OMA group metadata, and `last_modified` information for genomes.
- NCBI LinkOut resource generation script.

### Changed

- Domain annotations use the new Gene3D file format.
- Opening a database with a different major schema version raises an error; a minor difference only logs a warning. The schema version was bumped.
- HOG IDs always start with `HOG:`.
- Much faster computation of domain-covered sites.

### Fixed

- `resolve_oma_group`, values returned by the `Relations` model, and handling of missing data in standalone conversion.

## [0.5.0] - 2018-01-15

### Added

- Synteny conservation scores for homoeologs are stored in the database.
- More genome metadata (common name, URL, ...) and additional fields on the `Genome` model.
- Verification state of xrefs is stored in a new table column.
- Index on HOG levels.

### Changed

- Xref search uses a suffix array over lowercase xrefs.
- The `.idx` file is linked by a path relative to the main HDF5 file, and a missing `.idx` file is handled.

## [0.4.5] - 2018-01-11

### Changed

- Alternative splicing variants that share the same xrefs now resolve to the main variant.

## [0.4.4] - 2017-11-03

### Added

- `Database.main_isoforms()` to retrieve the main isoforms of a genome.
- Domain coverage statistics in the HDF5 database.

### Changed

- `biopython` added as a dependency.

### Fixed

- Python 2 compatibility fixes.

## [0.4.3] - 2017-10-06

### Added

- `Database.iter_members_of_hog_id()`.

### Fixed

- Bug in `member_of_fam`, and backwards compatibility of `ProteinModel.from_entry_nr`.

## [0.4.2] - 2017-10-04

### Added

- `parse_hog_id()` and `hog_members_from_hog_id()`.

### Fixed

- `search_xref` returns a proper int.

## [0.4.1] - 2017-10-02

### Changed

- The `ProteinEntry` model loads its entry lazily.

## [0.4.0] - 2017-09-26

### Added

- Sequence search (`SequenceSearch`) using a suffix array and k-mer index stored in the `.idx` file; `exact_search` has a `full_length` option.
- Gene Ontology and EC annotations in the database, with a `GeneOntology` lazy property and a GO FastMapper for function projection.
- More xref types, prefix ID matching, and xrefs stored in the suffix array/k-mer index.
- HOG domain prevalence, HOG metadata (size, longest sequence, top level), OMA group metadata, Fingerprints/Keywords, and Gene3D domains.
- New model properties: `Genome.EntryOff`, `Genome.totEntries`, `Genome.chromosomes`, `ProteinEntry.chromosome`, `locus_start`/end/strand, `md5` and GC content. `models.PairwiseRelation` was also added.
- Resolve genomes by NCBI taxon ID, and retrieve OMA group members by numeric group ID.
- `OrthoXMLSplitter`, Python computation of orthology relation types, and a basic database consistency-check framework.
- Summary statistics of the dataset.

### Changed

- Database schema 3.0; Newick serialization uses `[` `]` instead of `(` `)` in names.
- pyoma no longer depends on a specific HOG ID schema.
- Python 2 and 3 compatible (tested with tox).
- The `h5py` dependency was removed; `familyanalyzer` and other build dependencies moved to the `create_db` extra.

### Fixed

- Overflow bug in sequence offsets.
- Empty `orth.txt.gz` files during import.

## [0.3.1] - 2016-08-12

### Added

- Index on OMA groups.

### Fixed

- Conversion now includes the correct HOG content.
- Duplicates in the taxonomy table.
- Python 2 incompatibility in the input context manager.

## [0.3.0] - 2016-06-09

### Added

- Import of data from OMA standalone.
- Sanity checks to track plausibility between releases.
- Extended tests, with test data from git annex.

### Fixed

- Various small bugs.

## [0.2.1] - 2015-12-23

### Added

- Initial documentation of the library.
- Code to generate the HOG visualization data.

## [0.2.0] - 2015-12-18

- First tagged release.

[0.15.0]: https://github.com/DessimozLab/pyoma/compare/v0.14.2...v0.15.0
[0.14.2]: https://github.com/DessimozLab/pyoma/compare/v0.14.1...v0.14.2
[0.14.1]: https://github.com/DessimozLab/pyoma/compare/v0.14.0...v0.14.1
[0.14.0]: https://github.com/DessimozLab/pyoma/compare/v0.13.5...v0.14.0
[0.13.5]: https://github.com/DessimozLab/pyoma/compare/v0.13.4...v0.13.5
[0.13.4]: https://github.com/DessimozLab/pyoma/compare/v0.13.3...v0.13.4
[0.13.3]: https://github.com/DessimozLab/pyoma/compare/v0.13.2...v0.13.3
[0.13.2]: https://github.com/DessimozLab/pyoma/compare/v0.13.1...v0.13.2
[0.13.1]: https://github.com/DessimozLab/pyoma/compare/v0.13.0...v0.13.1
[0.13.0]: https://github.com/DessimozLab/pyoma/compare/v0.12.1...v0.13.0
[0.12.1]: https://github.com/DessimozLab/pyoma/compare/v0.12.0...v0.12.1
[0.12.0]: https://github.com/DessimozLab/pyoma/compare/v0.11.1...v0.12.0
[0.11.1]: https://github.com/DessimozLab/pyoma/compare/v0.11.0...v0.11.1
[0.11.0]: https://github.com/DessimozLab/pyoma/compare/v0.10.1...v0.11.0
[0.10.1]: https://github.com/DessimozLab/pyoma/compare/v0.10.0...v0.10.1
[0.10.0]: https://github.com/DessimozLab/pyoma/compare/0.9.4...v0.10.0
[0.9.4]: https://github.com/DessimozLab/pyoma/compare/0.9.3...0.9.4
[0.9.3]: https://github.com/DessimozLab/pyoma/compare/0.9.2...0.9.3
[0.9.2]: https://github.com/DessimozLab/pyoma/compare/0.9.1...0.9.2
[0.9.1]: https://github.com/DessimozLab/pyoma/compare/0.9.0...0.9.1
[0.9.0]: https://github.com/DessimozLab/pyoma/compare/0.8.0...0.9.0
[0.8.0]: https://github.com/DessimozLab/pyoma/compare/0.7.2...0.8.0
[0.7.2]: https://github.com/DessimozLab/pyoma/compare/0.7.1...0.7.2
[0.7.1]: https://github.com/DessimozLab/pyoma/compare/0.7.0...0.7.1
[0.7.0]: https://github.com/DessimozLab/pyoma/compare/0.6.4...0.7.0
[0.6.4]: https://github.com/DessimozLab/pyoma/compare/0.6.3...0.6.4
[0.6.3]: https://github.com/DessimozLab/pyoma/compare/0.6.2...0.6.3
[0.6.2]: https://github.com/DessimozLab/pyoma/compare/0.6.1...0.6.2
[0.6.1]: https://github.com/DessimozLab/pyoma/compare/0.6.0...0.6.1
[0.6.0]: https://github.com/DessimozLab/pyoma/compare/0.5.0...0.6.0
[0.5.0]: https://github.com/DessimozLab/pyoma/compare/0.4.5...0.5.0
[0.4.5]: https://github.com/DessimozLab/pyoma/compare/0.4.4...0.4.5
[0.4.4]: https://github.com/DessimozLab/pyoma/compare/0.4.3...0.4.4
[0.4.3]: https://github.com/DessimozLab/pyoma/compare/0.4.2...0.4.3
[0.4.2]: https://github.com/DessimozLab/pyoma/compare/0.4.1...0.4.2
[0.4.1]: https://github.com/DessimozLab/pyoma/compare/0.4.0...0.4.1
[0.4.0]: https://github.com/DessimozLab/pyoma/compare/0.3.1...0.4.0
[0.3.1]: https://github.com/DessimozLab/pyoma/compare/0.3.0...0.3.1
[0.3.0]: https://github.com/DessimozLab/pyoma/compare/0.2.1...0.3.0
[0.2.1]: https://github.com/DessimozLab/pyoma/compare/0.2.0...0.2.1
[0.2.0]: https://github.com/DessimozLab/pyoma/releases/tag/0.2.0
