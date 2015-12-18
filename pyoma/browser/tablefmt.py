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


class OrthoXmlHogTable(tables.IsDescription):
    Fam = tables.UInt32Col(pos=0)
    HogBufferOffset = tables.UInt32Col(pos=1)
    HogBufferLength = tables.UInt32Col(pos=2)


class ProteinTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    SeqBufferOffset = tables.UInt32Col(pos=2)
    SeqBufferLength = tables.UInt32Col(pos=3)
    OmaGroup = tables.UInt32Col(pos=4, dflt=0)
    OmaHOG = tables.StringCol(255, pos=5, dflt=b"")
    Chromosome = tables.StringCol(255, pos=6)
    LocusStart = tables.UInt32Col(pos=7)
    LocusEnd = tables.UInt32Col(pos=8)
    LocusStrand = tables.Int8Col(pos=9, dflt=1)
    AltSpliceVariant = tables.Int32Col(pos=10, dflt=0)
    CanonicalId = tables.StringCol(20, pos=11, dflt=b"")
    CDNABufferOffset = tables.UInt32Col(pos=12)
    CDNABufferLength = tables.UInt32Col(pos=13)
    MD5ProteinHash = tables.StringCol(32, pos=14)
    DescriptionOffset = tables.UInt32Col(pos=15)
    DescriptionLength = tables.UInt16Col(pos=16)
    SubGenome = tables.StringCol(1, pos=17, dflt=b"")


class LocusTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    Start = tables.UInt32Col(pos=2)
    End = tables.UInt32Col(pos=3)
    Strand = tables.Int8Col(pos=4)


class PairwiseRelationTable(tables.IsDescription):
    EntryNr1 = tables.UInt32Col(pos=0)
    EntryNr2 = tables.UInt32Col(pos=1)
    RelType = tables.EnumCol(
        tables.Enum({'1:1': 0, '1:n': 1, 'm:1': 2, 'm:n': 3,
                     'close paralog': 4, 'homeolog': 5, 'n/a': 6}),
        'n/a', base='uint8', pos=2)
    Score = tables.Float32Col(pos=3, dflt=-1)
    Distance = tables.Float32Col(pos=4, dflt=-1)
    AlignmentOverlap = tables.Float16Col(pos=5, dflt=-1)
    SyntenyConservationLocal = tables.Float16Col(pos=6, dflt=-1)
    SyntenyConservationChromosome = tables.Float16Col(pos=7, dflt=-1)


class XRefTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    XRefSource = tables.EnumCol(
        tables.Enum(['UniProtKB/SwissProt', 'UniProtKB/TrEMBL', 'n/a', 'EMBL',
                     'Ensembl Gene', 'Ensembl Transcript', 'Ensembl Protein',
                     'RefSeq', 'EntrezGene', 'GI', 'WikiGene', 'IPI',
                     'SourceID', 'SourceAC', 'PMP', 'NCBI', 'FlyBase']),
        'n/a', base='uint8', pos=2)
    XRefId = tables.StringCol(50, pos=3)


class GeneOntologyTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=1)
    TermNr = tables.UInt32Col(pos=2)
    Evidence = tables.StringCol(3, pos=3)
    Reference = tables.StringCol(255, pos=4)


class GenomeTable(tables.IsDescription):
    NCBITaxonId = tables.UInt32Col(pos=0)
    UniProtSpeciesCode = tables.StringCol(5, pos=1)
    TotEntries = tables.UInt32Col(pos=2)
    TotAA = tables.UInt32Col(pos=3)
    EntryOff = tables.UInt32Col(pos=4)
    SciName = tables.StringCol(255, pos=5)
    Release = tables.StringCol(255, pos=6)
    IsPolyploid = tables.BoolCol(pos=7)


class TaxonomyTable(tables.IsDescription):
    NCBITaxonId = tables.UInt32Col(pos=0)
    ParentTaxonId = tables.UInt32Col(pos=1)
    Name = tables.StringCol(255, pos=2)


class DomainTable(tables.IsDescription):
    EntryNr = tables.UInt32Col(pos=0)
    DomainId = tables.StringCol(20, pos=1)
    Coords = tables.StringCol(255, pos=2)


class DomainDescriptionTable(tables.IsDescription):
    DomainId = tables.StringCol(20, pos=0)
    Source = tables.StringCol(11, pos=1)
    Description = tables.StringCol(150, pos=2)
