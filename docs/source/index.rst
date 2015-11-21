.. GOparser documentation master file, created by
   sphinx-quickstart on Fri Sep 11 15:46:50 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GOparser |version|
==================

GOparser is a Python framework for working with Gene Ontology (GO) terms and
annotations. GOparser is free and open-source software (see
`License <license>`).

Main Features
-------------

- Efficient parsing of GO term information contained in ``go-basic.obo``
  from the
  `Gene Ontology Consortium <http://geneontology.org/page/download-ontology>`_.

- Efficient parsing of GO annotations contained in GAF files from EBI's
  `GO Annotation Database (UniProt-GOA) <http://www.ebi.ac.uk/GOA>`_.

- Filtering for annotations of protein-coding genes (using `genometools`).

- Easy and fast retrieval of all genes annotated with a particular GO term,
  or all GO terms a particular gene is annotated with.

- Annotations are fully propagated based on ``is_a`` and ``part_of``
  relations between GO terms.

- Cross-species support.

- Support for filtering annotations based on
  `evidence code <http://geneontology.org/page/guide-go-evidence-codes>`_.

Demo Notebooks
--------------

- `Demo.ipynb`__ (:download:`download <notebooks/Demo.ipynb>`): GOparser demo notebook

__ demo_notebook_

.. _demo_notebook: http://nbviewer.ipython.org/url/goparser.readthedocs.org/en/latest/_downloads/Demo.ipynb

.. This only links to the "latest" (master) version...no way to automatically switch to the "develop" version

Installation
------------

GOparser can be installed from `PyPI <https://pypi.python.org/pypi>`_ using
`pip <https://pip.pypa.io/en/stable/>`_:

.. code-block:: bash

    $ pip install goparser

Missing Features
----------------

- Visualizations (e.g., to show relationships between GO terms).
- Support for other GO relations, e.g., ``regulates`` and ``has_part``.

.. toctree::
    :maxdepth: 2
    :hidden:

    self
    Modules <modules>
    License <license>

.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`search`
.. * :ref:`modindex`
