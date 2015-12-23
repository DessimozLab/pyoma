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


Module documentation
####################

Here comes now the module doc:

.. automodule:: pyoma.browser.convert
    :members:

