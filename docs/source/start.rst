Getting Started
===============

.. _installation:

**Installation**
----------------
Python dependencies (defined in the `requirements <https://github.com/phac-nml/cladeomatic/blob/main/requirements.txt>`_ file, should be automatically installed when using conda or pip).

In addition to the python dependencies, Clade-o-Matic requires `Jellyfish 2.3.0 <https://github.com/gmarcais/Jellyfish/>`_, `psutil 5.9.1 <https://github.com/giampaolo/psutil>`_, `scikit-learn 1.1.1 <https://scikit-learn.org>`_, and `snp-dists 0.8.2 <https://github.com/tseemann/snp-dists>`_.
Install the latest released version from conda:

.. code-block:: console

    conda create -c bioconda -c conda-forge -n cladeomatic cladeomatic

Install using pip:

.. code-block:: console

    pip install cladeomatic

Install the latest master branch version directly from Github:

.. code-block:: console

    conda install jellyfish=2.30
    conda install psutil=5.9.1
    conda install scikit-learn=1.1.1
    conda install snp-dist=0.8.2
    pip install git+https://github.com/phac-nml/cladeomatic.git

.. _quickstart:

**Quick Start**
---------------
**Basic Usage**

If you run ``cladeomatic``, you should see the following usage statement:

.. code-block:: console

    Usage: cladeomatic <command> [options] <required arguments>

    To get minimal usage for a command use:
    cladeomatic command

    To get full help for a command use one of:
    cladeomatic command -h
    cladeomatic command --help

    Available commands:

    create     Identify population structure and develop typing scheme
    genotype   Call genotypes from a VCF and a scheme file
    benchmark  Test developed scheme using labeled samples and scheme*
    namer      Rename genotypes within a scheme*

    *in development


For further reference on the options available, refer the the :doc:`usage` guide.

**Create Scheme**
Option 1 - De novo tree-based

This mode will discover clades and lineages which meet membership size and SNP requirements.

Input requirements are:

* newick formatted tree
* VCF
* Reference sequence (.fasta / .gbk)
* Name of outgroup sequence
* Metadata file

``cladeomatic create --in_nwk tree.nwk  --in_var variants.vcf --in_meta metadata.txt --outdir scheme/ --root_name ref --reference ref.gbk``

Option 2 - Predefined groups

This mode will attempt to define a scheme based on a group manifest which meet membership size and SNP requirements.

Input requirements are:

* TSV formatted group file (sample_id, genotype)
* VCF
* Reference sequence (.fasta / .gbk)
* Name of outgroup sequence
* Metadata file

``cladeomatic create --in_groups groups.tsv --in_var variants.vcf --in_meta metadata.txt --outdir scheme/ --root_name ref --reference ref.gbk``


Outputs:

OutputFolderName

* {prefix}-altseq.fasta
* {prefix}-biohansel.fasta
* {prefix}-biohansel.meta.txt
* {prefix}-clades.info.txt
* {prefix}-dist.mat.txt
* {prefix}-extracted.kmers.txt
* {prefix}-filtered.vcf
* {prefix}-genotypes.distance.txt
* {prefix}-genotypes.raw.txt
* {prefix}-genotypes.selected.txt
* {prefix}-genotypes.supported.txt
* {prefix}-kmer.scheme.txt
* {prefix}-params.log
* {prefix}-sample.distances.html
* {prefix}-snps.scheme.txt
* {prefix}-snps.info.txt
* pseudo.seqs.fasta

