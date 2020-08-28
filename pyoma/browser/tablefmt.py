import tables

"""This module contains the definitions of the database tables
used in the browser database. Some of these tables are used
multiple times, e.g. the PairwiseRelationTable is used
for each genome pair.

From these table definitions one can easily extract the numpy
dtype that can hold the data:

   >>>tables.dtype_from_descr(HOGsTable)
   dtype([('Fam', '<i4'), ('ID', 'S255'), ('Level', 'S255')])
"""


class HOGsTable(tables.IsDescription):
    Fam = tables.Int32Col(pos=1)
    ID = tables.StringCol(255, pos=2)
    Level = tables.StringCol(255, pos=3)
    CompletenessScore = tables.Float32Col(pos=4, dflt=-1)
    ImpliedLosses = tables.Int32Col(pos=5, dflt=-1)
    NrMemberGenes = tables.Int32Col(pos=6, dflt=-1)
    IsRoot = tables.BoolCol(pos=7, dflt=False)


class OrthoXmlHogTable(tables.IsDescription):
    Fam = tables.UInt32Col(pos=0)
    HogBufferOffset = tables.Int64Col(pos=1)
    HogBufferLength = tables.UInt32Col(pos=2)
    HogAugmentedBufferOffset = tables.Int64Col(pos=3)
    HogAugmentedBufferLength = tables.UInt32Col(pos=4)


class AncestralSyntenyRels(tables.IsDescription):
    HogRow1 = tables.UInt32Col(pos=0)
    HogRow2 = tables.UInt32Col(pos=1)
    Weight = tables.Float16Col(pos=2)


class ProteinTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    SeqBufferOffset = tables.UInt64Col(pos=2)
    SeqBufferLength = tables.UInt32Col(pos=3)
    OmaGroup = tables.UInt32Col(pos=4, dflt=0)
    OmaHOG = tables.StringCol(255, pos=5, dflt=b"")
    Chromosome = tables.StringCol(255, pos=6)
    LocusStart = tables.UInt32Col(pos=7)
    LocusEnd = tables.UInt32Col(pos=8)
    LocusStrand = tables.Int8Col(pos=9, dflt=1)
    AltSpliceVariant = tables.Int32Col(pos=10, dflt=0)
    CanonicalId = tables.StringCol(20, pos=11, dflt=b"")
    CDNABufferOffset = tables.UInt64Col(pos=12)
    CDNABufferLength = tables.UInt32Col(pos=13)
    MD5ProteinHash = tables.StringCol(32, pos=14)
    DescriptionOffset = tables.UInt32Col(pos=15)
    DescriptionLength = tables.UInt16Col(pos=16)
    SubGenome = tables.StringCol(1, pos=17, dflt=b"")
    RootHogUpstream = tables.Int32Col(pos=18, dflt=-1)
    RootHogDownStream = tables.Int32Col(pos=19, dflt=-1)


class ProteinCacheInfo(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=0)
    NrPairwiseOrthologs = tables.UInt32Col(pos=1)
    NrHogInducedPWOrthologs = tables.UInt32Col(pos=2)
    NrHogInducedPWParalogs = tables.UInt32Col(pos=3)
    NrOMAGroupOrthologs = tables.UInt32Col(pos=4)
    NrAnyOrthologs = tables.UInt32Col(pos=5)


class LocusTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    Start = tables.UInt32Col(pos=2)
    End = tables.UInt32Col(pos=3)
    Strand = tables.Int8Col(pos=4)


class PairwiseRelationTable(tables.IsDescription):
    EntryNr1 = tables.UInt32Col(pos=0)
    EntryNr2 = tables.UInt32Col(pos=1)
    RelType = tables.EnumCol(
        tables.Enum(
            {
                "1:1": 0,
                "1:n": 1,
                "m:1": 2,
                "m:n": 3,
                "close paralog": 4,
                "homeolog": 5,
                "n/a": 6,
            }
        ),
        "n/a",
        base="uint8",
        pos=2,
    )
    Score = tables.Float32Col(pos=3, dflt=-1)
    Distance = tables.Float32Col(pos=4, dflt=-1)
    AlignmentOverlap = tables.Float16Col(pos=5, dflt=-1)
    SyntenyConservationLocal = tables.Float16Col(pos=6, dflt=-1)
    Confidence = tables.Float16Col(pos=7, dflt=-1)


class XRefTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    XRefSource = tables.EnumCol(
        tables.Enum(
            {
                "UniProtKB/SwissProt": 0,
                "UniProtKB/TrEMBL": 10,
                "Ensembl Protein": 20,
                "Ensembl Gene": 25,
                "Ensembl Transcript": 30,
                "RefSeq": 40,
                "EntrezGene": 50,
                "FlyBase": 60,
                "WormBase": 65,
                "EnsemblGenomes": 70,
                "NCBI": 75,
                "EMBL": 80,
                "SourceID": 95,
                "SourceAC": 100,
                "HGNC": 105,
                "Gene Name": 110,
                "Synonym": 115,
                "Protein Name": 120,
                "ORF Name": 125,
                "Ordered Locus Name": 130,
                "PDB": 148,
                "Swiss Model": 150,
                "STRING": 151,
                "neXtProt": 152,
                "Bgee": 153,
                "EPD": 154,
                "ChEMBL": 156,
                "GlyConnect": 157,
                "SwissPalm": 158,
                "DisGeNET": 159,
                "WikiGene": 160,
                "IPI": 240,
                "GI": 241,
                "n/a": 255,
            }
        ),  # last line: deprecated systems
        "n/a",
        base="uint8",
        pos=2,
    )
    XRefId = tables.StringCol(50, pos=3)
    Verification = tables.EnumCol(
        tables.Enum({"exact": 0, "unchecked": 2, "modified": 4}),
        "unchecked",
        base="uint8",
        pos=4,
    )


class GeneOntologyTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    TermNr = tables.UInt32Col(pos=2)
    Evidence = tables.StringCol(3, pos=3)
    Reference = tables.StringCol(255, pos=4)


class ECTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    ECacc = tables.StringCol(16, pos=2)


class GenomeTable(tables.IsDescription):
    NCBITaxonId = tables.Int32Col(pos=0)
    UniProtSpeciesCode = tables.StringCol(5, pos=1)
    TotEntries = tables.UInt32Col(pos=2)
    TotAA = tables.UInt32Col(pos=3)
    EntryOff = tables.UInt32Col(pos=4)
    SciName = tables.StringCol(255, pos=5)
    CommonName = tables.StringCol(64, pos=6)
    SynName = tables.StringCol(64, pos=7)
    Release = tables.StringCol(128, pos=8)
    Url = tables.StringCol(255, pos=9)
    Source = tables.StringCol(255, pos=10)
    Date = tables.Time32Col(pos=11)
    IsPolyploid = tables.BoolCol(pos=12)


class TaxonomyTable(tables.IsDescription):
    NCBITaxonId = tables.Int32Col(pos=0)
    ParentTaxonId = tables.Int32Col(pos=1)
    Name = tables.StringCol(255, pos=2)


class DomainTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=0)
    DomainId = tables.StringCol(20, pos=1)
    Coords = tables.StringCol(255, pos=2)


class DomainDescriptionTable(tables.IsDescription):
    DomainId = tables.StringCol(20, pos=0)
    Source = tables.StringCol(11, pos=1)
    Description = tables.StringCol(150, pos=2)


class HOGDomainArchPrevalenceTable(tables.IsDescription):
    Fam = tables.UInt32Col(pos=0)
    ReprEntryNr = tables.UInt32Col(pos=1)
    ReprEntryLength = tables.UInt32Col(pos=2)
    TopLevel = tables.StringCol(255, pos=3)
    FamSize = tables.UInt32Col(pos=4)
    PrevCount = tables.UInt32Col(pos=5)


class HOGDomainPresenceTable(tables.IsDescription):
    Offset = tables.UInt32Col(pos=0)
    DomainId = tables.StringCol(20, pos=1)


class GroupsizeHistogram(tables.IsDescription):
    Size = tables.UInt32Col(pos=0)
    Count = tables.UInt32Col(pos=1)


class OmaGroupTable(tables.IsDescription):
    GroupNr = tables.UInt32Col(pos=0)
    Fingerprint = tables.StringCol(7, pos=1)
    KeywordOffset = tables.UInt32Col(pos=2)
    KeywordLength = tables.UInt16Col(pos=3)
    NrMembers = tables.UInt16Col(pos=4)


class RootHOGMetaTable(tables.IsDescription):
    FamNr = tables.UInt32Col(pos=0)
    KeywordOffset = tables.UInt32Col(pos=1)
    KeywordLength = tables.UInt32Col(pos=2)
    FamDataJsonOffset = tables.UInt64Col(pos=3)
    FamDataJsonLength = tables.UInt32Col(pos=4)
