pyOMA - python library to interact with OMA hdf5 files
======================================================

TODO: extend readme.

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
