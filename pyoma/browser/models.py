from __future__ import division

import collections
import numpy
import time


def format_sciname(sci, short=False):
    p = set(
        [
            sci.find(x)
            for x in [
                "(",
                "serogroup",
                "serotype",
                "serovar",
                "biotype",
                "subsp",
                "pv.",
                "bv.",
            ]
        ]
    )
    if sci.startswith("Escherichia coli"):
        p.add(sci.find("O"))
    p.discard(-1)
    p = min(p) if len(p) > 0 else len(sci)
    return {"species": sci[0:p], "strain": sci[p:]}


class LazyProperty(object):
    """Decorator to evaluate a property only on access.

    Compute the attribute value and caches it in the instance.
    Python Cookbook (Denis Otkidach) http://stackoverflow.com/users/168352/denis-otkidach
    This decorator allows you to create a property which can be computed once and
    accessed many times."""

    def __init__(self, method, name=None):
        # record the unbound-method and the name
        self.method = method
        self.name = name or method.__name__
        self.__doc__ = method.__doc__

    def __get__(self, inst, cls):
        if inst is None:
            return self
        # compute, cache and return the instance's attribute value
        result = self.method(inst)
        # setattr redefines the instance's attribute so this doesn't get called again
        setattr(inst, self.name, result)
        return result


class KeyWrapper(object):
    """
    Enables the use of functions, e.g. bisect, with a key function.
    """

    def __init__(self, it, key):
        self.it = it
        self.key = key

    def __getitem__(self, i):
        return self.key(self.it[i])

    def __len__(self):
        return len(self.it)


class Singleton(type):
    """A meta-class to enforce a Singleton, e.g. a class that can be
    instantiated only exactly once.

    Modified from Python Cookbook, 3rd Edition, p 357ff.

    :Example:

        class Foo(metaclass=Singleton):
            def __init__(self):
                pass  #This part is executed only once
    """

    def __init__(self, *args, **kwargs):
        self.__instance = None
        super(Singleton, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        if self.__instance is None:
            self.__instance = super(Singleton, self).__call__(*args, **kwargs)
        return self.__instance


class ProteinEntry(object):
    """Model for a protein object

    This class provides an easy to use interface for a given protein
    form the database.

    If instantiated with an entry_nr only, no data is loaded until a
    property or method is accessed. Properties that need to access
    additional data or loaded lazily and are cached in the object
    (but not kept after deletion of object)."""

    def __init__(self, db, e):
        self._stored_entry = e
        self._db = db

    @LazyProperty
    def _entry(self):
        return (
            self._db.entry_by_entry_nr(self._stored_entry)
            if isinstance(self._stored_entry, (int, numpy.integer))
            else self._stored_entry
        )

    @classmethod
    def from_entry_nr(cls, db, eNr):
        # e = db.entry_by_entry_nr(eNr)
        return cls(db, int(eNr))

    @property
    def entry_nr(self):
        return int(self._entry["EntryNr"])

    @property
    def locus_start(self):
        return int(self._entry["LocusStart"])

    @property
    def locus_end(self):
        return int(self._entry["LocusEnd"])

    @property
    def strand(self):
        return int(self._entry["LocusStrand"])

    @LazyProperty
    def exons(self):
        return ExonStructure.from_entry_nr(self._db, self.entry_nr)

    @property
    def nr_exons(self):
        return int(len(self.exons))

    @property
    def oma_group(self):
        return int(self._entry["OmaGroup"])

    @property
    def oma_hog(self):
        return self._entry["OmaHOG"].decode()

    @property
    def chromosome(self):
        return self._entry["Chromosome"].decode()

    @property
    def canonicalid(self):
        return self._entry["CanonicalId"].decode()

    @LazyProperty
    def xrefs(self):
        return self._db.id_mapper["Xref"].map_entry_nr(self._entry["EntryNr"])

    @property
    def sequence_md5(self):
        return self._entry["MD5ProteinHash"].decode()

    @LazyProperty
    def genome(self):
        g = self._db.id_mapper["OMA"].genome_of_entry_nr(self._entry["EntryNr"])
        return Genome(self._db, g)

    @LazyProperty
    def omaid(self):
        return self._db.id_mapper["OMA"].map_entry_nr(self._entry["EntryNr"])

    @LazyProperty
    def cdna(self):
        return self._db.get_cdna(self._entry).decode()

    @property
    def gc_content(self):
        cdna = self.cdna
        cnts = list(map(cdna.count, "GCAT"))
        try:
            return sum(cnts[0:2]) / sum(cnts)
        except ZeroDivisionError:
            return 0

    @LazyProperty
    def sequence(self):
        return self._db.get_sequence(self._entry).decode()

    @property
    def sequence_length(self):
        return int(self._entry["SeqBufferLength"]) - 1

    @LazyProperty
    def description(self):
        return self._db.get_description(self._entry).decode()

    @LazyProperty
    def ec_numbers(self):
        return self._db.get_ec_annotations(self.entry_nr)

    @property
    def subgenome(self):
        return self._entry["SubGenome"].decode()

    @LazyProperty
    def hog_family_nr(self):
        from .db import Singleton as HOGSingleton

        try:
            fam = self._db.hog_family(self._entry)
        except HOGSingleton:
            fam = 0
        return fam

    @property
    def is_main_isoform(self):
        return bool(
            self._entry["AltSpliceVariant"] == 0
            or self._entry["AltSpliceVariant"] == self._entry["EntryNr"]
        )

    @LazyProperty
    def alternative_isoforms(self):
        return [
            ProteinEntry(self._db, e)
            for e in self._db.get_splicing_variants(self._entry)
            if e["EntryNr"] != self.entry_nr
        ]

    def get_main_isoform(self):
        if self.is_main_isoform:
            return self
        else:
            return ProteinEntry(self._db, self._entry["AltSpliceVariant"])

    def __repr__(self):
        return "<{}({}, {})>".format(self.__class__.__name__, self.entry_nr, self.omaid)

    def __len__(self):
        return self.sequence_length


class Genome(object):
    """Model of a genome/proteome

    This model provides information about a genome. It is instantiated with
    row of the the /Genome table from the hdf5 file."""

    def __init__(self, db, g):
        self._genome = g
        self._db = db

    @property
    def ncbi_taxon_id(self):
        """returns the ncbi taxonomy id of the genome"""
        return int(self._genome["NCBITaxonId"])

    @property
    def uniprot_species_code(self):
        """returns the uniprot mnemonic species code"""
        return self._genome["UniProtSpeciesCode"].decode()

    @property
    def sciname(self):
        """returns the scientific name of the genome"""
        return self._genome["SciName"].decode()

    @property
    def common_name(self):
        """returns the common name of the genome. If not set,
        the empty string is returned"""
        try:
            return self._genome["CommonName"].decode()
        except ValueError:
            return ""

    @property
    def synonym_name(self):
        """returns the synonym name"""
        return self._genome["SynName"].decode()

    @LazyProperty
    def species_and_strain_as_dict(self):
        return format_sciname(self.sciname)

    def species(self):
        return self.species_and_strain_as_dict["species"]

    def strain(self):
        return self.species_and_strain_as_dict["strain"]

    @property
    def url(self):
        return self._genome["Url"].decode()

    @property
    def source(self):
        return self._genome["Source"].decode()

    @property
    def release(self):
        return self._genome["Release"].decode()

    @property
    def last_modfied_timestamp(self):
        return self._genome["Date"]

    @property
    def last_modified(self):
        return self.modification_date("%Y-%b-%d")

    def modification_date(self, fmt):
        if self._db.db_schema_version >= (3, 2):
            return time.strftime(fmt, time.localtime(self.last_modfied_timestamp))
        else:
            return "n/a"

    @property
    def nr_entries(self):
        """returns the number of protein entries of the genome in the database"""
        return int(self._genome["TotEntries"])

    @LazyProperty
    def nr_genes(self):
        """returns the number of genes of the genome in the database"""
        return self._db.count_main_isoforms(self.uniprot_species_code)

    @property
    def entry_nr_offset(self):
        return int(self._genome["EntryOff"])

    @LazyProperty
    def kingdom(self):
        # TODO: store directly in db
        return self._db.tax.get_parent_taxa(self._genome["NCBITaxonId"])[-1][
            "Name"
        ].decode()

    @property
    def is_polyploid(self):
        """returns whether or not this is a polyploid genome"""
        return self._genome["IsPolyploid"]

    @LazyProperty
    def lineage(self):
        """returns the lineage of the genome."""
        return [
            lev["Name"].decode()
            for lev in self._db.tax.get_parent_taxa(self._genome["NCBITaxonId"])
        ]

    @LazyProperty
    def chromosomes(self):
        chrs = collections.defaultdict(list)
        entry_tab = self._db.get_hdf5_handle().get_node("/Protein/Entries")
        for row in entry_tab.where(
            "(EntryNr > {}) & (EntryNr <= {})".format(
                self.entry_nr_offset, self.entry_nr_offset + self.nr_entries
            )
        ):
            chrs[row["Chromosome"].decode()].append(row["EntryNr"])
        return chrs

    def approx_chromosome_length(self, chromosome):
        """method to retrieve the approximate length of a given chromosome.
        The value corresponds to the end locus coordinate of the last gene
        on the chromosome.

        :param chromosome: the chromosome of interest
        :type chromosome: str, bytes
        :raises ValueError: if chromosome does not exist for this genome"""
        if isinstance(chromosome, bytes):
            chr = chromosome
        elif isinstance(chromosome, str):
            chr = chromosome.encode("utf-8")
        else:
            chr = "{}".format(chromosome).encode("utf-8")
        query = "(EntryNr > {}) & (EntryNr <= {}) & (Chromosome == {!r})".format(
            self.entry_nr_offset, self.entry_nr_offset + self.nr_entries, chr
        )
        tab = self._db.get_hdf5_handle().get_node("/Protein/Entries")
        chr_len = int(max((row["LocusEnd"] for row in tab.where(query))))
        return chr_len

    def __repr__(self):
        return "<{}({}, {})>".format(
            self.__class__.__name__, self.uniprot_species_code, self.ncbi_taxon_id
        )

    def __len__(self):
        return self.nr_entries


class OmaGroup(object):
    """OmaGroup object model

    The OmaGroup model can be instantiated with a group nr or a fingerprint.
    The (meta-)data for the group will be loaded lazily if properties are accessed."""

    def __init__(self, db, og):
        self._stored_group = og
        self._db = db

    @LazyProperty
    def _group(self):
        if isinstance(self._stored_group, dict) and "group_nr" in self._stored_group:
            return self._stored_group
        else:
            return self._db.oma_group_metadata(
                self._db.resolve_oma_group(self._stored_group)
            )

    @property
    def group_nbr(self):
        """numeric representation of the OmaGroup"""
        return int(self._group["group_nr"])

    @property
    def fingerprint(self):
        """fingerprint of the OmaGroup"""
        return self._group["fingerprint"]

    @property
    def keyword(self):
        """inferred keyword of the OmaGroup"""
        return self._group["keywords"]

    @LazyProperty
    def members(self):
        """returns a list of member proteins as :class:`ProteinEntry` objects"""
        return [
            ProteinEntry(self._db, pe)
            for pe in self._db.oma_group_members(self.group_nbr)
        ]

    @property
    def nr_member_genes(self):
        """number of genes/proteins belonging to he group."""
        size = self._group["size"]
        if size < 0:
            # fallback in case it is not yet stored in db
            size = len(self.members)
        return size

    def __len__(self):
        return self.nr_member_genes

    def __repr__(self):
        return "<{}({}, {})>".format(
            self.__class__.__name__, self.group_nbr, self.fingerprint
        )


class HOG(object):
    """HOG object model

    This object stores information about a HOG at a certain level and provides
    convenince methods / properties to load various additional data related to
    the HOG from the underlying database.

    The model can be instantiated with a *hog* argument, that can be either a
    roothog family interger numer, or a valid Hog-ID (as string/byte). Optionally,
    the *level* argument can be used to select the sub-hog at the given taxonomic
    range. Up on access of any property / method of the
    created object, the model will lazily load the
    :class:`pyoma.browser.tablefmt.HogLevel` array for the specified hog / level
    arguments.

    :param db: the underlying database object
    :type db: :class:`pyoma.browser.db.Database`
    :param hog: the hog id or the HogLevel numpy array instance
    :type hog: int, str, bytes, numpy.void
    :param level: desired level of (sub-)HOG
    :type level: str, bytes, optional"""

    def __init__(self, db, hog, level=None):
        self._stored_hog = hog
        self._stored_level = level
        self._db = db

    @LazyProperty
    def _hog(self):
        if isinstance(self._stored_hog, (int, numpy.integer)):
            # load hog at root
            return self._db.get_hog(
                self._db.format_hogid(self._stored_hog), level=self._stored_level
            )
        elif isinstance(self._stored_hog, (str, bytes)):
            # load hog using ID
            return self._db.get_hog(self._stored_hog, self._stored_level)
        else:
            return self._stored_hog

    @property
    def fam(self):
        """roothog id as integer"""
        return int(self._hog["Fam"])

    @property
    def hog_id(self):
        """HogID of the HOG"""
        return self._hog["ID"].decode()

    @property
    def level(self):
        """taxonomic range of the HOG"""
        return self._hog["Level"].decode()

    @property
    def nr_member_genes(self):
        """number of extent genes belonging to the HOG"""
        return int(self._hog["NrMemberGenes"])

    @property
    def is_root(self):
        """whether or not this is the deepest taxonomic range for
        this (sub-)HOG."""
        return bool(self._hog["IsRoot"])

    @property
    def completeness_score(self):
        """returns the completness score of the HOG. It corresponds to the fraction of
        genomes that have at least one gene in this HOG."""
        return float(self._hog["CompletenessScore"])

    def __len__(self):
        return self.nr_member_genes

    def __repr__(self):
        return "<{}({}, {})>".format(self.__class__.__name__, self.hog_id, self.level)

    @LazyProperty
    def keyword(self):
        """returns the inferred keywords for the HOG. It is based on the descriptions
        and IDs of all the member genes"""
        return self._db.get_roothog_keywords(self.fam)

    @LazyProperty
    def members(self):
        """returns a list of :class:`pyoma.browser.models.ProteinEntry` instances with all
        the proteins / genes belonging to this HOG."""
        return [
            ProteinEntry(self._db, pe)
            for pe in self._db.member_of_hog_id(self.hog_id, level=self.level)
        ]

    @LazyProperty
    def parent_hogs(self):
        """returns a list of the parent HOGs."""
        return self._db.get_parent_hogs(self.hog_id, self.level)


class PairwiseRelation(object):
    def __init__(self, db, relation):
        self._relation = relation
        self._db = db

    @property
    def distance(self):
        return float(self._relation["Distance"])

    @property
    def score(self):
        return float(self._relation["Score"])

    @property
    def alignment_overlap(self):
        return float(self._relation["AlignmentOverlap"])

    @property
    def synteny_conservation_local(self):
        return float(self._relation["SyntenyConservationLocal"])

    @property
    def confidence(self):
        return float(self._relation["Confidence"])

    @LazyProperty
    def rel_type(self):
        if not isinstance(self._relation["RelType"], str):
            type_map = self._db._get_pw_tab(
                self._relation["EntryNr1"], "VPairs"
            ).get_enum("RelType")
            return type_map(self._relation["RelType"])
        else:
            return self._relation["RelType"]

    @LazyProperty
    def entry_1(self):
        return ProteinEntry(
            self._db, self._db.entry_by_entry_nr(self._relation["EntryNr1"])
        )

    @LazyProperty
    def entry_2(self):
        return ProteinEntry(
            self._db, self._db.entry_by_entry_nr(self._relation["EntryNr2"])
        )


class GeneOntologyAnnotation(object):
    def __init__(self, db, anno):
        self.db = db
        self.anno = anno

    @LazyProperty
    def term(self):
        return self.db.gene_ontology.term_by_id(self.anno["TermNr"])

    @property
    def evidence(self):
        return self.anno["Evidence"].decode()

    @property
    def reference(self):
        return self.anno["Reference"].decode()

    @property
    def entry_nr(self):
        return int(self.anno["EntryNr"])

    @LazyProperty
    def aspect(self):
        from .geneontology import GOAspect

        return GOAspect.to_string(self.term.aspect)


class ExonStructure(object):
    def __init__(self, db, exons):
        self._stored = exons
        self._db = db

    @LazyProperty
    def _exons(self):
        return (
            self._db.get_exons(self._stored)
            if isinstance(self._stored, int)
            else self._stored
        )

    @classmethod
    def from_entry_nr(cls, db, eNr):
        return cls(db, int(eNr))

    def _iter_exons(self):
        if len(self._exons) > 0 and self._exons["Strand"][0] < 0:
            self._exons[::-1].sort(order="Start")
        else:
            self._exons.sort(order="Start")
        for exon in self._exons:
            yield Exon(exon)

    def __len__(self):
        return len(self._exons)

    def __repr__(self):
        return "<{}(entry_nr={}, nr_exons={})>".format(
            self.__class__.__name__, self._exons[0]["EntryNr"], len(self)
        )

    def __str__(self):
        exs = list(str(e) for e in self._iter_exons())
        if len(exs) > 1:
            return "join({})".format(", ".join(exs))
        elif len(exs) == 1:
            return exs[0]
        else:
            return "n/a"

    def as_list_of_dict(self):
        return [
            {"start": e.start, "end": e.end, "strand": e.strand}
            for e in self._iter_exons()
        ]


class Exon(object):
    def __init__(self, exon):
        self.exon = exon

    @property
    def start(self):
        return int(self.exon["Start"])

    @property
    def end(self):
        return int(self.exon["End"])

    @property
    def strand(self):
        return "+" if self.exon["Strand"] > 0 else "-"

    def __str__(self):
        if self.exon["Strand"] < 0:
            template = "complement({}..{})"
        else:
            template = "{}..{}"
        return template.format(self.exon["Start"], self.exon["End"])
