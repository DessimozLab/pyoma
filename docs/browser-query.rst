Querying the OMA Browser
========================

The module :mod:`pyoma.browser.db` is the main entrypoint to work with data from 
an OMA browser instance in python. The first thing to do is to get an instance of 
the  database object:

.. code-block:: python

    import pyoma.browser.db
    db = pyoma.browser.db.Database('./OmaServer.h5')

which will open the hdf5 file you specify. This object provides access to lot's of data,
especially the orthology predictions, gene neighborhood, sequences and meta-data:

.. autoclass:: pyoma.browser.db.Database
    :members:

Further, the db object keeps references to the idmapper dispatcher, an id resolver and the 
taxonomy object through the attributes `id_mapper`, `id_resolver` and `tax`.

Resolving an ID
###############

The id resolver can be used to convert an crossreference of a protein id into the numeric 
entry_id of the protein. The class provides the following interface:

.. autoclass:: pyoma.browser.db.IDResolver
    :members:

The `resolve` method can be used to map any xref to the entry number, also the 
OMA ids (e.g. HUMAN02123).


Accessing the species phylogeny
###############################

OMA provides access to the species phylogeny it used to infer the HOGs. Currently this is the 
NCBI Taxonomy. The :class:`Taxonomy` has the follwing interface that can be used to query 
the hirarchy:

.. autoclass:: pyoma.browser.db.Taxonomy
    :members:

Essentially this boils down to be able to get the parental lineages of a query genome.


Mapping from and to Crossreferences
###################################

To convert an internal entry_nr to any standard identifier one needs the id_mapper functionality.
The attribute :attr:`id_mapper` in the database object points to an factory that depending on the 
requested type, returns the associated mapper object.

The following mappers exist so far:

  *OMA* : mapping from and to OMA ids (e.g. HUMAN13233). Internally also needed to determine 
          the genome of an protein entry.

  *Xref* : mapping from and to any external ids. 

  *UniProt* : subset of crossreferences that point to UniProt or SwissProt entries.

  *Linkout* : subset of crossreference that point to UniProt, Ensembl and Entrez. 
              The method :meth:`xreftab_to_dict` will create an additional
              key *url* that points to the entry's website of the respecive resource.


The individual types can be used as accessors of the facotry to get the respecive mapper 
object, e.g.

.. code-block:: python

    db.id_mapper['OMA']
    db.id_mapper['XRef']


The OMA mapper provides the following interface:

.. autoclass:: pyoma.browser.db.OmaIdMapper
    :members:

Similarly for the XRef id mapper object:

.. autoclass:: pyoma.browser.db.XrefIdMapper
    :members:


Working with GeneOntology annotations
#####################################

Gene Annotations for a specific protein (or for a
range of entries) in the database can be loaded with the
:meth:`pyoma.browser.db.Database.get_gene_ontology_annotations`.

In order to traverse the ontology itself (e.g. to get also the
parent terms of a specific annotation, a
:class:`pyoma.browser.gene_ontology.GeneOntology` object needs
to be created. The easiest way for this is by using the
property :meth:`pyoma.browser.db.Database.gene_ontology` of the
database.



.. doctest::

    >>> import pyoma.browser.db
    >>> db = pyoma.browser.db.Database('tests/browser/TestDb.h5')
    >>> annos = db.get_gene_ontology_annotations(11522)
    >>> annos
    array([(11522,  3735, b'ISS', b'PMID:12368867'),
           (11522,  3735, b'IEA', b'GO_REF:002'),
           (11522,  5622, b'IEA', b'GO_REF:002'),
           (11522,  5840, b'IEA', b'GO_REF:002'),
           (11522,  5840, b'IEA', b'GO_REF:038'),
           (11522,  6412, b'ISS', b'PMID:12368867'),
           (11522,  6412, b'IEA', b'GO_REF:002'),
           (11522, 22627, b'ISS', b'PMID:12368867'),
           (11522, 30529, b'IEA', b'GO_REF:038')],
          dtype=[('EntryNr', '<u4'), ('TermNr', '<u4'), ('Evidence', 'S3'), ('Reference', 'S255')])
    >>> term = db.gene_ontology.term_by_id(annos[0]['TermNr'])
    >>> print("{0} - {0.name}".format(term))
    GO:0003735 - structural constituent of ribosome

Working with GO can also be simplified by using working with the
:class:`pyoma.browser.model.GeneOntology` model:

.. code-block:: python

    >>> import pyoma.browser.models
    >>> annos = [pyoma.browser.models.GeneOntologyAnnotation(db, anno) for anno in annos]
    >>> annos[0].evidence
    'ISS'




