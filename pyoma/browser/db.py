from builtins import chr
from builtins import range
from builtins import object
from builtins import zip
import itertools

import tables
import numpy
import numpy.lib.recfunctions
import re
import json
import os
import collections
import logging

logger = logging.getLogger(__name__)


def count_elements(iterable):
    """return the number of elements in an iterator in the most efficient way.

    Be aware that for unbound iterators, this method won't terminate!
    :param iterable: an iterable object.
    """
    counter = itertools.count()
    collections.deque(zip(iterable, counter), maxlen=0)  # (consume at C speed)
    return next(counter)


class Database(object):
    """This is the main interface to the oma database. Queries
    will typically be issued by methods of this object. Typically
    the result of queries will be :py:class:`numpy.recarray` objects."""
    EXPECTED_DB_SCHEMA = "2.0"

    def __init__(self, db):
        if isinstance(db, str):
            self.db = tables.open_file(db, 'r')
        elif isinstance(db, tables.File):
            self.db = db
        else:
            raise ValueError(str(db) + ' is not a valid database type')

        try:
            db_version = self.db.get_node_attr('/', 'db_schema_version')
        except AttributeError:
            db_version = "1.0"

        if db_version != self.EXPECTED_DB_SCHEMA:
            raise DBVersionError('Unsupported database version: {} != {} ({})'
                                 .format(db_version, self.EXPECTED_DB_SCHEMA, self.db.filename))

        self.id_resolver = IDResolver(self)
        self.id_mapper = IdMapperFactory(self)
        self.tax = Taxonomy(self.db.root.Taxonomy.read())

        self._re_fam = re.compile(r'HOG:(?P<fam>\d{7,})')

    def get_hdf5_handle(self):
        """return the handle to the database hdf5 file"""
        return self.db

    def ensure_entry(self, entry):
        """This method allows to use an entry or an entry_nr.

        If necessary it will load the entry from the entry_nr,
        otherwise returning the same object again.

        :param entry: the entry_nr of a protein to be loaded or a
            protein entry."""
        try:
            t = entry['AltSpliceVariant']
            return entry
        except (TypeError, AttributeError, IndexError):
            if isinstance(entry, (int, numpy.number)):
                return self.entry_by_entry_nr(entry)
            raise TypeError('Invalid type to retrieve an Entry')
        except Exception:
            raise TypeError('Invalid type to retrieve an Entry')

    def entry_by_entry_nr(self, entry_nr):
        """Returns the entry from the /Protein/Entries table
        corresponding to entry_nr.

        :param int entry_nr: a numeric identifier for the protein
            entry"""
        entry = self.db.root.Protein.Entries[entry_nr - 1]
        if entry['EntryNr'] != entry_nr:
            logger.warning('EntryNr {} not at position {}. Using index instead'.format(entry_nr, entry_nr - 1))
            entry = self.db.root.Protein.Entries.read_where(
                'EntryNr == {:d}'.format(entry_nr))
            if len(entry) != 1:
                raise ValueError("there are {} entries with entry_nr {}".format(len(entry), entry_nr))
            entry = entry[0]
        return entry

    def _get_vptab(self, entry_nr):
        return self._get_pw_tab(entry_nr, 'VPairs')

    def _get_pw_tab(self, entry_nr, subtab):
        genome = self.id_mapper['OMA'].genome_of_entry_nr(entry_nr)['UniProtSpeciesCode'].decode()
        return self.db.get_node('/PairwiseRelation/{}/{}'.format(genome, subtab))

    def count_vpairs(self, entry_nr):
        vptab = self._get_vptab(entry_nr)
        try:
            cnt = count_elements(vptab.where('(EntryNr1=={:d})'.format(entry_nr)))
        except (TypeError, ValueError):
            cnt = 0
        return cnt

    def _get_pw_data(self, entry_nr, tab):
        dat = tab.read_where('(EntryNr1=={:d})'.format(entry_nr))
        typ = tab.get_enum('RelType')
        res = numpy.lib.recfunctions.append_fields(
            dat[['EntryNr1', 'EntryNr2', 'Score', 'Distance']],
            names='RelType',
            data=[typ(x) for x in dat['RelType']],
            usemask=False)
        return res

    def get_vpairs(self, entry_nr):
        """returns the verified pairs of a query protein.

        This method returns an instance of a :class:`numpy.recarray` class
        containing the verified pairs of a query protein entry.
        The returned array contains columns with EntryNr1 and EntryNr2 to
        identify the pair together with RelType (indicating the subtype of
        orthology), the alignment score and the distance. The score and
        distance will be set to -1 if unknown.

        :param int entry_nr: the numeric entry_nr of the query protein."""
        vp_tab = self._get_vptab(entry_nr)
        return self._get_pw_data(entry_nr, vp_tab)

    def get_within_species_paralogs(self, entry_nr):
        """returns the within species paralogs of a given entry

        This method returns a :class:`numpy.recarray` instance
        containing the close paralogs. Close paralogs are within
        species paralogs that are inparalogs to at least one
        ortholog of the query gene in OMA.

        The returned array contains columns with EntryNr1 and EntryNr2 to
        identify the pair together with RelType (indicating the subtype of
        paralogy), the alignment score and the distance. The score and
        distance will be set to -1 if unknown.

        :param int entry_nr: the numeric entry_id of the query protein"""
        within_species_paralogs = self._get_pw_tab(entry_nr, 'within')
        return self._get_pw_data(entry_nr, within_species_paralogs)

    def neighbour_genes(self, entry_nr, window=1):
        """Returns neighbor genes around a query gene.

        This method returns a tuple containing a numpy recarray with
        gene entries located around the query gene, and an index
        pointing to the query gene. The genes are sorted according to
        their position on the chromosome.

        The *windows* parameter specifies the number of genes up- and
        downstream of the query gene that should be reported. Note
        that the actual number can be smaller if the query gene is close
        to a chromosome start or end.

        :param entry_nr: the entry number of the query gene
        :param window: the number of neighboring genes on each
                           side to return"""
        if window <= 0 or not isinstance(window, int):
            raise ValueError('windows parameters must be a positive integer value')

        dat = self.entry_by_entry_nr(entry_nr)
        target_chr = dat['Chromosome']
        genome_range = self.id_mapper['OMA'].genome_range(entry_nr)
        f = 5
        data = self.db.root.Protein.Entries.read_where(
            '(EntryNr >= {:d}) & (EntryNr <= {:d}) & '
            '(Chromosome == {!r}) & '
            '((AltSpliceVariant == 0) |'
            ' (AltSpliceVariant == EntryNr))'.format(
                max(genome_range[0], entry_nr - f * window),
                min(genome_range[1], entry_nr + f * window),
                target_chr))
        data.sort(order=['EntryNr'])
        idx = data['EntryNr'].searchsorted(entry_nr)
        res = data[max(0, idx - window):min(len(data), idx + window + 1)]
        idx = res['EntryNr'].searchsorted(entry_nr)
        return res, idx

    def hog_family(self, entry):
        entry = self.ensure_entry(entry)
        m = self._re_fam.match(entry['OmaHOG'].decode())
        if m is None:
            raise Singleton(entry)
        return int(m.group('fam'))

    def hog_levels_of_fam(self, fam_nr):
        """get all taxonomic levels covered by a family.

        The family coresponds to the toplevel numeric id of a HOG,
        i.e. for HOG:002421 the fam_nr should be 2421. If a HOG
        covers a certain level more than once, it will be returned
        several times.

        :param fam_nr: the numeric id of the family (== Toplevel HOG)
        """
        return self.db.root.HogLevel.read_where(
            '(Fam=={})'.format(fam_nr))['Level']

    def get_subhogids_at_level(self, fam_nr, level):
        """get all the hog ids within a given family at a given taxonomic
        level of interest.

        After a duplication in an ancestor lineage, there exists multiple
        sub-hogs for any taxonomic level after the duplication. This method
        allows to get the list of hogids at the requested taxonomic level.

        E.g. assume in family 1 (HOG:0000001) there has been a duplication
        between Eukaryota and Metazoa. this method would return for
        get_subhogids_at_level(1, 'Eukaryota') --> ['HOG:0000001']
        and for
        get_subhogids_at_level(1, 'Metazoa') --> ['HOG:0000001.1a', 'HOG:0000001.1b']

        :param fam_nr: the numeric family id
        :param level: the taxonomic level of interest"""
        lev = level if isinstance(level, bytes) else level.encode('ascii')
        return self.db.root.HogLevel.read_where(
            '(Fam=={}) & (Level=={!r})'.format(fam_nr, lev))['ID']

    def member_of_hog_id(self, hog_id, level=None):
        """return an array of protein entries which belong to a given hog_id.

        E.g. if hog_id = 'HOG122.1a', the method returns all the proteins that
        have either exactly this hog id or an inparalogous id such a HOG122.1a.4b.2a

        If you are only interested in the members of a specific lineage (identified
        through its taxonomic range), you can pass the taxonomic range as an
        additional argument. Only the proteins of genomes belonging to this clade
        will be returned. Otherwise, all proteins with having this specific hog_id
        will be returned.

        :param str hog_id: the requested hog_id.
        :param level: the taxonomic level of interest
        :type level: str or None

        :return: a numpy.array with the protein entries belonging to the requested hog.
        :rtype: :class:`numpy.ndarray`

        :Note: Even if you obtained a certain hog_id using
               :py:meth:`get_subhogids_at_level`
               using a certain level, if you do not specify the level in
               :meth:`member_of_hog_id` again, you will likely get proteins from other
               clades. Only if it happens that the deepest level of the hog_id
               coincides with the taxonomic range of interest, the two will be identical.
        """
        hog_range = self._hog_lex_range(hog_id)
        # get the proteins which have that HOG number
        memb = self.db.root.Protein.Entries.read_where(
            '({!r} <= OmaHOG) & (OmaHOG < {!r})'.format(*hog_range))
        if level is not None:
            memb = [x for x in memb if level.encode('ascii') in self.tax.get_parent_taxa(
                self.id_mapper['OMA'].genome_of_entry_nr(x['EntryNr'])['NCBITaxonId'])['Name']]

        return memb

    def member_of_fam(self, fam):
        """returns an array of protein entries which belong to a given fam"""
        if not isinstance(fam, (int, numpy.number)):
            raise ValueError('expect a numeric family id')
        return self.member_of_hog_id('HOG:{:07d}'.format(fam))

    def hog_members(self, entry, level):
        """get hog members with respect to a given taxonomic level.

        The method will return a list of protein entries that are all
        member of the same hog with respect to the taxonomic range
        of interest.

        :param entry: an entry or entry_nr of a query protein
        :param level: the taxonomic level of interest"""
        query = self.ensure_entry(entry)
        queryFam = self.hog_family(query)
        hoglev = None
        for hog_candidate in self.db.root.HogLevel.where(
                '(Fam == {:d}) & (Level == {!r})'.format(queryFam, level.encode('ascii'))):
            if query['OmaHOG'].startswith(hog_candidate['ID']):
                hoglev = hog_candidate
                break
        if hoglev is None:
            raise ValueError(u'Level "{0:s}" undefined for query gene'.format(level))
        # get the entries which have this hogid (or a sub-hog)
        members = self.member_of_hog_id(hoglev['ID'])
        if level != 'LUCA':
            # last, we need to filter the proteins to the tax range of interest
            members = [x for x in members if level.encode('ascii') in self.tax.get_parent_taxa(
                self.id_mapper['OMA'].genome_of_entry_nr(x['EntryNr'])['NCBITaxonId'])['Name']]
            if query not in members:
                raise ValueError(u"Level '{0:s}' undefined for query gene".format(level))
        return members

    def get_orthoxml(self, fam):
        """returns the orthoxml of a given toplevel HOG family

        :param fam: numeric id of requested toplevel hog"""
        idx = self.db.root.OrthoXML.Index.read_where('Fam == {:d}'.format(fam))
        if len(idx) < 1:
            raise ValueError('cannot retrieve orthoxml for {}'.format(fam))
        idx = idx[0]
        return self.db.root.OrthoXML.Buffer[
               idx['HogBufferOffset']:idx['HogBufferOffset'] + idx['HogBufferLength']].tostring()

    def _hog_lex_range(self, hog):
        """return the lexographic range of a hog.

        This can be used to search of sub-hogs which are nested in
        the query hog. The semantics is such that
        _hog_lex_range[0] <= hog < _hog_lex_range[1].
        This is equivalent to say that a sub-hog starts with the
        query hog."""
        hog_str = hog.decode() if isinstance(hog, bytes) else hog
        return hog_str.encode('ascii'), (hog_str[0:-1] + chr(1 + ord(hog_str[-1]))).encode('ascii')

    def oma_group_members(self, group_nr):
        """get the member entries of an oma group.

        This method returns a numpy array of protein entries that form
        an oma group. If the numeric group id is invalid (not positive
        integer value), an `InvalidId` Exception is raised.

        :param int group_nr: numeric oma group id"""
        if not isinstance(group_nr, int) or group_nr <= 0:
            raise InvalidId('Invalid id of oma group')
        members = self.db.root.Protein.Entries.read_where('OmaGroup=={:d}'.format(group_nr))
        return members

    def get_sequence(self, entry):
        """get the protein sequence of a given entry as a string

        :param entry: the entry or entry_nr for which the sequence is requested"""
        entry = self.ensure_entry(entry)
        seqArr = self.db.get_node('/Protein/SequenceBuffer')
        seq = seqArr[entry['SeqBufferOffset']:entry['SeqBufferOffset'] + entry['SeqBufferLength'] - 1]
        return seq.tostring()

    def get_cdna(self, entry):
        """get the protein sequence of a given entry as a string"""
        entry = self.ensure_entry(entry)
        seqArr = self.db.get_node('/Protein/CDNABuffer')
        seq = seqArr[entry['CDNABufferOffset']:entry['CDNABufferOffset'] + entry['CDNABufferLength'] - 1]
        return seq.tostring()

    def get_description(self, entry):
        entry = self.ensure_entry(entry)
        descArr = self.db.get_node('/Protein/DescriptionBuffer')
        desc = descArr[entry['DescriptionOffset']:entry['DescriptionOffset'] + entry['DescriptionLength']]
        return desc.tostring()

    def get_release_name(self):
        return str(self.db.get_node_attr('/', 'oma_version'))

    def get_domains(self, entry_nr):
        try:
            return self.db.root.Annotations.Domains.read_where('EntryNr == {:d}'.format(entry_nr))
        except ValueError as e:
            raise InvalidId('require a numeric entry id, got {}'.format(entry_nr))


class OmaIdMapper(object):
    def __init__(self, db):
        self.genome_table = db.get_hdf5_handle().root.Genome.read()
        self._entry_off_keys = self.genome_table.argsort(order=('EntryOff'))
        self._genome_keys = self.genome_table.argsort(
            order=('UniProtSpeciesCode'))
        self._omaid_re = re.compile(r'(?P<genome>[A-Z][A-Z0-9]{4})(?P<nr>\d+)')
        self._db = db

    def genome_of_entry_nr(self, e_nr):
        """returns the genome code belonging to a given entry_nr"""
        idx = self.genome_table['EntryOff'].searchsorted(
            e_nr - 1, side='right',
            sorter=self._entry_off_keys)
        return self.genome_table[self._entry_off_keys[idx - 1]]

    def map_entry_nr(self, entry_nr):
        genome = self.genome_of_entry_nr(entry_nr)
        return "{0:s}{1:05d}".format(genome['UniProtSpeciesCode'].decode(),
                                     entry_nr - genome['EntryOff'])

    def genome_from_UniProtCode(self, code):
        code = code.encode('ascii')
        idx = self.genome_table['UniProtSpeciesCode'].searchsorted(
            code, sorter=self._genome_keys)
        genome = self.genome_table[self._genome_keys[idx]]
        if genome['UniProtSpeciesCode'] != code:
            raise UnknownSpecies('{} is unknown'.format(code))
        return genome

    def omaid_to_entry_nr(self, omaid):
        """returns the internal numeric entrynr from a
        UniProtSpeciesCode+nr id. this is the inverse
        function of 'map_entry_nr'."""
        match = self._omaid_re.match(omaid)
        if match is None:
            raise InvalidOmaId(omaid)
        code, nr = match.group('genome'), int(match.group('nr'))
        genome = self.genome_from_UniProtCode(code)
        if nr <= 0 or nr > genome['TotEntries']:
            raise InvalidOmaId(omaid)
        return genome['EntryOff'] + int(match.group('nr'))

    def genome_range(self, query):
        """returns the internal range of EntryNr associated with
        'query'. 'query' can be either a numeric id of a protein
        or a UniProtSpeciesCode of a genome. If 'query' is unknown
        by the database, an InvalidOmaId exception is raised.

        The return range is a tuple of length two, and the numbers
        indicated the *inclusive* boundaries, e.g. (1,5) indicates
        that the entries 1,2,3,4 and 5 belong to the query species"""
        if isinstance(query, (int, numpy.integer)):
            genome_row = self.genome_of_entry_nr(query)
            if query <= 0 or query > genome_row['EntryOff'] + genome_row['TotEntries']:
                raise InvalidOmaId(query)
        else:
            genome_row = self.genome_from_UniProtCode(query)
        return (genome_row['EntryOff'] + 1,
                genome_row['EntryOff'] + genome_row['TotEntries'],)

    def species_ordering(self, root=None):
        """get ordering of the genomes with respect to taxonomy.

        This method returns a linear ordering of all the genomes with
        respect to their lineage, i.e. genomes that are evolutionary
        "close" to each other appear close in the ordering.
        Optionally, one can give a root genome, that will be the species
        the ordering is going to start with.

        :param root: UniProtSpeciesCode of the root genome.
        :returns: a list of species codes in the correct order."""
        if root is None:
            root = self.genome_table[0]['UniProtSpeciesCode']
        root_genome = self.genome_from_UniProtCode(root)
        lins = {g['UniProtSpeciesCode']: [lev['Name'] for lev in self._db.tax.get_parent_taxa(g['NCBITaxonId'])][::-1]
                for g in self.genome_table}
        root_lin = lins[root_genome['UniProtSpeciesCode']]
        sort_key = {}
        for g, lin_g in lins.items():
            for k in range(min(len(root_lin), len(lin_g))):
                if root_lin[k] != lin_g[k]:
                    k -= 1
                    break
            sort_key[g] = (-k, lin_g)
        sorted_genomes = sorted(list(sort_key.keys()), key=lambda g: sort_key[g])
        return {g.decode(): v for v, g in enumerate(sorted_genomes)}


class IDResolver(object):
    def __init__(self, db):
        entry_nr_col = db.get_hdf5_handle().root.Protein.Entries.cols.EntryNr
        self.max_entry_nr = entry_nr_col[int(entry_nr_col.index[-1])]
        self._db = db

    def _from_numeric(self, e_id):
        nr = int(e_id)
        if not 0 < nr <= self.max_entry_nr:
            raise InvalidId('{0:d} out of protein range: {1:s}'.format(nr, e_id))
        return nr

    def _from_omaid(self, e_id):
        return self._db.id_mapper['OMA'].omaid_to_entry_nr(e_id)

    def search_xrefs(self, e_id):
        """search for all xrefs. TODO: what happens if xref is ambiguous?"""
        res = set([x['EntryNr'] for x in self._db.id_mapper['XRef'].search_xref(e_id)])
        if len(res) == 0:
            raise InvalidId(e_id)
        elif len(res) > 1:
            raise NotImplemented('Cross-ref "{}" is ambiguous'.format(e_id))
        return int(res.pop())

    def resolve(self, e_id):
        """maps an id to the entry_nr of the current OMA release."""
        try:
            nr = self._from_numeric(e_id)
        except ValueError:
            try:
                nr = self._from_omaid(e_id)
            except (InvalidOmaId, UnknownSpecies) as e:
                nr = self.search_xrefs(e_id)
        return nr


class Taxonomy(object):
    """Taxonomy provides an interface to navigate the taxonomy data.

    The input data is the same as what is stored in the Database in
    table "/Taxonomy"."""

    def __init__(self, data):
        if not isinstance(data, numpy.ndarray):
            raise ValueError('Taxonomy expects a numpy table.')
        self.tax_table = data
        self.taxid_key = self.tax_table.argsort(order=('NCBITaxonId'))
        self.parent_key = self.tax_table.argsort(order=('ParentTaxonId'))
        self._load_valid_taxlevels()

    def _load_valid_taxlevels(self):
        forbidden_chars = re.compile(r'[^A-Za-z. -]')
        try:
            with open(os.environ['DARWIN_BROWSERDATA_PATH'] + '/TaxLevels.drw') as f:
                taxStr = f.read()
            tax_json = json.loads(("[" + taxStr[14:-3] + "]").replace("'", '"'))
            self.all_hog_levels = frozenset([t.encode('ascii') for t in
                                             tax_json if forbidden_chars.search(t) is None])
        except (IOError, KeyError):
            self.all_hog_levels = frozenset([l for l in self.tax_table['Name']
                                             if forbidden_chars.search(l.decode()) is None])

    def _table_idx_from_numeric(self, tid):
        i = self.tax_table['NCBITaxonId'].searchsorted(
            tid, sorter=self.taxid_key)
        idx = self.taxid_key[i]
        if self.tax_table[idx]['NCBITaxonId'] != tid:
            raise InvalidTaxonId(u"{0:d} is an invalid/unknown taxonomy id".format(tid))
        return idx

    def _get_root_taxon(self):
        i1 = self.tax_table['ParentTaxonId'].searchsorted(0, sorter=self.parent_key)
        i2 = self.tax_table['ParentTaxonId'].searchsorted(0, sorter=self.parent_key, side='right')
        if i2 - i1 == 0:
            raise DBConsistencyError('Not a single root in Taxonomy: {}'
                                     .format(self.tax_table[self.parent_key[i1]]))
        elif i2 - i1 == 1:
            res = self.tax_table[self.parent_key[i1]]
        else:
            res = numpy.array([(0, -1, b'LUCA')], dtype=self.tax_table.dtype)[0]
        return res

    def _taxon_from_numeric(self, tid):
        idx = self._table_idx_from_numeric(tid)
        return self.tax_table[idx]

    def _direct_children_taxa(self, tid):
        i = self.tax_table['ParentTaxonId'].searchsorted(tid, sorter=self.parent_key)
        idx = []
        while i < len(self.parent_key) and self.tax_table[self.parent_key[i]]['ParentTaxonId'] == tid:
            idx.append(self.parent_key[i])
            i += 1
        return self.tax_table.take(idx)

    def get_parent_taxa(self, query):
        """Get array of taxonomy entries leading towards the
        root of the taxonomy.

        :param query: the starting taxonomy level"""
        idx = []
        parent = query
        count = 0
        while parent != 0:
            i = self._table_idx_from_numeric(parent)
            idx.append(i)
            tmp = self.tax_table[i]['ParentTaxonId']
            if tmp == parent:
                raise InvalidTaxonId(u"{0:d} has itself as parent".format(tmp))
            parent = tmp
            count += 1
            if count > 100:
                raise InvalidTaxonId(u"{0:d} exceeds max depth of 100. Infinite recursion?".format(query))
        return self.tax_table.take(idx)

    def _get_taxids_from_any(self, it, skip_missing=True):
        if not isinstance(it, numpy.ndarray):
            try:
                it = numpy.fromiter(it, dtype='i4')
            except ValueError:
                it = numpy.fromiter(it, dtype='S255')
        if it.dtype.type is numpy.string_:
            try:
                ns = self.name_key
            except AttributeError:
                ns = self.name_key = self.tax_table.argsort(order='Name')
            idxs = self.tax_table['Name'].searchsorted(it, sorter=ns)
            idxs = numpy.clip(idxs, 0, len(ns) - 1)
            taxs = self.tax_table[ns[idxs]]
            keep = taxs['Name'] == it
            if not skip_missing and not keep.all():
                raise KeyError('not all taxonomy names could be found')
            res = taxs['NCBITaxonId'][keep]
        else:
            res = it
        return res

    def get_induced_taxonomy(self, members, collapse=True):
        """Extract the taxonomy induced by a given set of `members`.

        This method allows to extract the part which is induced bay a
        given set of levels and leaves that should be part of the
        new taxonomy. `members` must be an iterable, the levels
        must be either numeric taxids or scientific names.

        :param iter members: an iterable containing the levels
            and leaves that should remain in the new taxonomy. can be
            either axonomic ids or scientific names.
        :param bool collapse: whether or not levels with only one child
            should be skipped or not. This defaults to True"""
        taxids_to_keep = numpy.sort(self._get_taxids_from_any(members))
        idxs = numpy.searchsorted(self.tax_table['NCBITaxonId'], taxids_to_keep, sorter=self.taxid_key)
        idxs = numpy.clip(idxs, 0, len(self.taxid_key) - 1)
        subtaxdata = self.tax_table[self.taxid_key[idxs]]
        if not numpy.alltrue(subtaxdata['NCBITaxonId'] == taxids_to_keep):
            raise KeyError('not all levels in members exists in this taxonomy')

        updated_parent = numpy.zeros(len(subtaxdata), 'bool')
        for i, cur_tax in enumerate(taxids_to_keep):
            if updated_parent[i]:
                continue
            # get all the parents and check which ones we keep in the new taxonomy.
            parents = self.get_parent_taxa(cur_tax)['NCBITaxonId']
            mask = numpy.in1d(parents, taxids_to_keep)
            # find the position of them in subtaxdata (note: subtaxdata and
            # taxids_to_keep have the same ordering).
            new_idx = taxids_to_keep.searchsorted(parents[mask])
            taxids = taxids_to_keep[new_idx]
            # parent taxid are ncbitaxonids shifted by one position!
            parents = numpy.roll(taxids, -1)
            parents[-1] = 0
            subtaxdata['ParentTaxonId'][new_idx] = parents
            updated_parent[new_idx] = True

        if collapse:
            nr_children = collections.defaultdict(int)
            for p in subtaxdata['ParentTaxonId']:
                nr_children[p] += 1
            rem = [p for (p, cnt) in nr_children.items() if cnt == 1 and p != 0]
            if len(rem) > 0:
                idx = taxids_to_keep.searchsorted(rem)
                return self.get_induced_taxonomy(numpy.delete(taxids_to_keep, idx))
        return Taxonomy(subtaxdata)

    def newick(self):
        """Get a Newick representation of the Taxonomy

        Note: as many newick parsers do not support quoted labels,
        the method instead replaces spaces with underscores."""

        def _rec_newick(node):
            children = []
            for child in self._direct_children_taxa(node['NCBITaxonId']):
                children.append(_rec_newick(child))

            if len(children) == 0:
                return node['Name'].decode().replace(' ', '_')
            else:
                t = ",".join(children)
                return '(' + t + ')' + node['Name'].decode().replace(' ', '_')

        return _rec_newick(self._get_root_taxon())

    def as_dict(self):
        """Encode the Taxonomy as a nested dict.

         This representation can for example be used to serialize
         a Taxonomy in json format."""

        def _rec_phylogeny(node):
            res = {'name': node['Name'].decode(), 'id': int(node['NCBITaxonId'])}
            children = []
            for child in self._direct_children_taxa(node['NCBITaxonId']):
                children.append(_rec_phylogeny(child))
            if len(children) > 0:
                res['children'] = children
            return res

        return _rec_phylogeny(self._get_root_taxon())


class InvalidTaxonId(Exception):
    pass


class DBVersionError(Exception):
    pass


class DBConsistencyError(Exception):
    pass


class InvalidId(Exception):
    pass


class InvalidOmaId(InvalidId):
    pass


class UnknownIdType(Exception):
    pass


class UnknownSpecies(Exception):
    pass


class Singleton(Exception):
    def __init__(self, entry, msg=None):
        super(Singleton, self).__init__(msg)
        self.entry = entry


class IdMapperFactory(object):
    def __init__(self, db_obj):
        self.db = db_obj
        self.mappers = {}

    def __getitem__(self, idtype):
        return self.get_mapper(idtype)

    def get_mapper(self, idtype):
        try:
            mapper = self.mappers[idtype]
        except KeyError:
            try:
                mapper = globals()[str(idtype).title() + 'IdMapper'](self.db)
                self.mappers[idtype] = mapper
            except KeyError:
                raise UnknownIdType('{} is unknown'.format(str(idtype)))
        return mapper


class XrefIdMapper(object):
    def __init__(self, db):
        self._db = db
        self.xref_tab = db.get_hdf5_handle().get_node('/XRef')
        self.xrefEnum = self.xref_tab.get_enum('XRefSource')
        self.idtype = frozenset(list(self.xrefEnum._values.keys()))

    def map_entry_nr(self, entry_nr):
        """returns the set of XRef entries associated with the query
        protein.

        The types of XRefs that are returned depends on the idtype
        class member variable. In the base-class, idtype contains
        all valid xref types. Typically, subclasses of XrefIdMapper
        will change this set.

        :param entry_nr: the numeric id of the query protein."""
        res = set([row[:] for row in self.xref_tab.where('EntryNr=={:d}'.format(entry_nr))
                   if row['XRefSource'] in self.idtype])
        return res

    def _combine_query_values(self, field, values):
        parts = ['({}=={})'.format(field, z) for z in values]
        return '|'.join(parts)

    def map_many_entry_nrs(self, entry_nrs):
        """map several entry_nrs with as few db queries as possible
        to their cross-references. The function returns a
        :class:`numpy.recarray` containing all fields as defined in
        the table.

        :param entry_nrs: a list with numeric protein entry ids"""
        mapped_junks = []
        junk_size = 32 - len(self.idtype)  # respect max number of condition variables.
        source_condition = self._combine_query_values('XRefSource', self.idtype)
        for start in range(0, len(entry_nrs), junk_size):
            condition = "({}) & ({})".format(
                self._combine_query_values('EntryNr',
                                           entry_nrs[start:start + junk_size]),
                source_condition)
            mapped_junks.append(self.xref_tab.read_where(condition))
        return numpy.lib.recfunctions.stack_arrays(
            mapped_junks,
            usemask=False)

    def search_xref(self, xref):
        """identify proteins associcated with `xref`.

        The crossreferences are limited to the types in the class
        member `idtype`. In the base class, all types are valid
        xrefs. The method returns a :class:`numpy.recarry` defined
        for the XRef table with all entries pointing to `xref`.

        :param xref: an xref to be located"""
        res = self.xref_tab.read_where('XRefId=={!r}'.format(xref.encode('utf-8')))
        if len(res) > 0 and len(self.idtype) < len(self.xrefEnum):
            res = res[numpy.in1d(res['XRefSource'], list(self.idtype))]
        return res

    def xreftab_to_dict(self, tab):
        """convert a xreftable to a dictionary per entry_nr.

        All rows in `tab` are converted into a nested dictionary
        where the outer key is a protein entry number and the
        inner key the xref source type.

        :param tab: a :class:`numpy.recarray` corresponding to XRef
            table definition to be converted"""
        xrefdict = collections.defaultdict(dict)
        for row in tab:
            try:
                typ = self.xrefEnum._values[row['XRefSource']]
            except IndexError:
                logger.warn('invalid XRefSource value in {}'.format(row))
                continue
            if typ not in xrefdict[row['EntryNr']]:
                xrefdict[row['EntryNr']][typ] = {'id': row['XRefId']}
        return xrefdict


class UniProtIdMapper(XrefIdMapper):
    def __init__(self, db):
        super(UniProtIdMapper, self).__init__(db)
        self.idtype = frozenset([self.xrefEnum[z]
                                 for z in ['UniProtKB/SwissProt', 'UniProtKB/TrEMBL']])


class LinkoutIdMapper(XrefIdMapper):
    def __init__(self, db):
        super(LinkoutIdMapper, self).__init__(db)
        self.idtype = frozenset([self.xrefEnum[z]
                                 for z in ['UniProtKB/SwissProt', 'UniProtKB/TrEMBL',
                                           'Ensembl Protein', 'EntrezGene']])

    def url(self, typ, id_):
        # TODO: improve url generator in external module with all xrefs
        url = None
        id_ = id_.decode()
        if typ.startswith('UniProtKB'):
            url = 'http://uniprot.org/uniprot/{}'.format(id_)
        elif typ == 'EntrezGene':
            url = 'http://www.ncbi.nlm.nih.gov/gene/{}'.format(id_)
        elif typ.startswith('Ensembl'):
            url = 'http://ensembl.org/id/{}'.format(id_)
        return url

    def xreftab_to_dict(self, tab):
        xref = super(LinkoutIdMapper, self).xreftab_to_dict(tab)
        for d in list(xref.values()):
            for typ, elem in list(d.items()):
                elem['url'] = self.url(typ, elem['id'])
        return xref


class DomainNameIdMapper(object):
    def __init__(self, db):
        self.domain_src = db.get_hdf5_handle().root.Annotations.DomainDescription.read()
        self.domain_src.sort(order='DomainId')

    def _get_dominfo(self, domain_id):
        idx = self.domain_src['DomainId'].searchsorted(domain_id)
        if self.domain_src[idx]['DomainId'] != domain_id:
            raise KeyError("no domain info available for {}".format(domain_id))
        return self.domain_src[idx]

    def get_info_dict_from_domainid(self, domain_id):
        info = self._get_dominfo(domain_id)
        return {'name': info['Description'].decode(), 'source': info['Source'].decode(),
                'domainid': domain_id.decode()}
