.. cladeomatic documentation master file, created by
   sphinx-quickstart on Fri Apr 21 12:34:10 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Clade-o-matic's documentation!
=========================================
Identification of population structure based on genomic data is a routine problem in molecular epidemiology and numerous approaches have been developed to partition large collections of sequences into clusters at different levels of granularity and resolution based on the goals of the researchers and underlying genetic diversity of the dataset. There is a confusing variety of terms in use to describe these units of genetic relatedness depending on the level of resolution provided and the specific taxon: clades, clones, clades, groups, genotypes, haplotypes, lineages, and subtypes.

Clade-o-matic is a phylogenetic approach for identification of hierarchal genotypes based on canonical SNPs that are exclusive and conserved within phylogenetic clades. It works best with organisms that have low rates of recombination and limited genetic diversity, but it is possible to analyze species-wide phylogenies. A major distinction of Clade-o-matic from other population structure tools such as rheirbaps is that it provides the distinguishing features for each clade in the form of a SNP typing scheme with unique kmers for identification of these SNPs from raw sequencing data. Furthermore, Clade-o-matic does not use a fixed set of levels for hierarchy as there will be uneven groupings between different clades. It uses the underlying distance distributions to compress the hierarchy into the minimum set of levels that retain the tree structure. The granularity of the scheme can be changed by specifiying different numbers of SNPs required to support a clade along with the number of members assigned to a clade for it to be considered valid.

.. note::

   This project is under active development.

Check out the :doc:`start` section for the Quick Start guide, including how to
:ref:`install <installation>` and :ref:`run <quickstart>` Clade-o-matic.

.. toctree::
   start
   :maxdepth: 3

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
