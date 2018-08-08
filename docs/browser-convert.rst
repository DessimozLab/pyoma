Convert a darwin oma browser instance to python OmaServer.h5
=============================================================


.. note::
    This code has been extracted from the browser/pyoma repository and moved
    in the pyoma library.

The :mod:`pyoma.browser.convert` submodule is needed to convert the data from an darwin
instance into a python compatible hdf5 database. The steps to create a darwin
oma browser are explained in detail on the `Wiki <http://lab.dessimoz.org/wiki/oma_browser_release>`_

As of now, the documentation is very small. Essentially, you can convert the database
with a single command, e.g. 

.. code-block:: sh

    bin/importdata.py --release /path/to/release OmaServer.h5

which will create the OmaServer.h5 file in the release path. If the optional argument is not 
specified, the release is extracted from the environment variable *DARWIN_BROWSERDATA_PATH*.


Note on CATH Domains
####################
CATH domains are available for uniprot's reference proteomes. The pipeline
to compute CATH domains is available and can be used to also predict domains
for sequences not annotated so far. The required steps are listed here:

  - TODO: list steps.


NCBI Linkout files
##################

Since 2017 we OMA is an accepted Linkout resource on the NCBI websites. For each new release
the data needs to be converted and uploaded to NCBI's FTP server. Integration from there
happens automatically every week or so. Details can be found on the
`NCBI Linkout page <https://www.ncbi.nlm.nih.gov/projects/linkout/doc/nonbiblinkout.html>`_

The script `bin/ncbi_linkout` does the conversion and also uploads the
resulting files (each step can also be done separately, note the command
line options for that). The password and user can be found on the group's
intranet wiki page.


Module documentation
####################

Here comes now the module doc:

.. automodule:: pyoma.browser.convert
    :members:

Splitting an OrthoXML file
--------------------------

A useful tool might be OrthoXMLSplitter that can be used also outside of the conversion
procedure. It has the following interface.

.. autoclass:: pyoma.browser.OrthoXMLSplitter.OrthoXMLSplitter
    :members:
    :special-members: __call__

