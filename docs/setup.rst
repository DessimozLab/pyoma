Installing pyoma library
========================

Currently, the only option to install the pyoma library is from the
group's git-repository

Installing from the git server
------------------------------

You need access to the git repositories. You can check with the following command:

.. code-block:: sh

    ssh gitolite@lab.dessimoz.org -p 2222

A list of the repositories you have access to should be listed. If you don't have
access, please send your public ssh key to Christophe or Adrian and request access.

Now, clone the zoo library and then install it into your current python environment:

.. code-block:: sh

    git clone ssh://gitolite@lab.dessimoz.org:2222/pyoma
    
    # install with resolved dependencies
    pip install -r requirements
    # otherwise, you might also use
    pip install -e .
    

Afterwards, you can use the functions by importing them into your python session, e.g.

.. code-block:: python
    
    import pyoma.browser.db

which contains the code to query the hdf5 files of an OMA Browser instance

