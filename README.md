pyOMA - python library to interact with OMA hdf5 files
======================================================

PyOMA is a python library that enables access to the hdf5 database of the
[OMA project]. Report problems to the [github issue] tracker.

Installing pyoma
----------------

    pip install pyoma

for support to build oma hdf5 files from either oma standalone runs or production files use

    pip install pyoma[create_db]

Quick tour
----------


    import pyoma.browser.db
    db = pyoma.browser.db.Database("OmaServer.h5")
    db.get_release_name()

Documentation
-------------

The documentation is available on https://zoo.cs.ucl.ac.uk/doc/pyoma

[OMA project]: https://omabrowser.org
[github issue]: https://github.com/DessimozLab/pyoma/issues
