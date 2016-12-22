import tables


class GenomeFmt(tables.IsDescription):
    NCBITaxonId = tables.UInt32Col(pos=0)
    UniProtSpeciesCode = tables.StringCol(5, pos=1)
    TotEntries = tables.UInt32Col(pos=2)
    TotAA = tables.UInt32Col(pos=3)
    EntryOff = tables.UInt32Col(pos=4)
    SciName = tables.StringCol(255, pos=5)
    Release = tables.StringCol(255, pos=6)


class TaxonomyFmt(tables.IsDescription):
    NCBITaxonId = tables.UInt32Col(pos=0)
    ParentTaxonId = tables.UInt32Col(pos=1)
    Name = tables.StringCol(255, pos=2)


class BestMatchesFmt(tables.IsDescription):
    EntryNr1 = tables.UInt32Col(pos=0)
    EntryNr2 = tables.UInt32Col(pos=1)
    Score = tables.Float32Col(pos=2)
    PamDistance = tables.Float32Col(pos=3)
    PamVariance = tables.Float32Col(pos=4)
    SumLengths = tables.UInt16Col(pos=5)
    isSP = tables.BoolCol(pos=6, dflt=True)
    isVP = tables.BoolCol(pos=7, dflt=False)
    isCanonicalSplicing = tables.BoolCol(pos=8, dflt=True)
