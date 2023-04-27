Getting Started
===============

.. _installation:

**Installation**
----------------
Python dependencies (defined in the `requirements <https://github.com/phac-nml/cladeomatic/blob/main/requirements.txt>`_ file, should be automatically installed when using conda or pip).

In addition to the python dependencies, Clade-o-Matic requires `Jellyfish 2.3.0 <https://github.com/gmarcais/Jellyfish/>`_ and `snp-dists 0.8.2 <https://github.com/tseemann/snp-dists>`_.
Install the latest released version from conda:

.. code-block:: console

    conda create -c bioconda -c conda-forge -n cladeomatic cladeomatic

Install using pip:

.. code-block:: console

    pip install cladeomatic

Install the latest master branch version directly from Github:

.. code-block:: console

    conda install jellyfish snp-dists
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
    benchmark  Test developed scheme using labeled samples and scheme
    namer      Rename genotypes within a scheme


For further reference on the options available, refer the the :doc:`usage` guide.

**Create Scheme**

Option 1 - De novo tree-based

This mode will discover clades and lineages which meet membership size and SNP requirements.

Input requirements are:

* newick formatted tree
* Reference (Outgroup) sequence (.fasta / .gbk)
* VCF (Must use the same reference sequence as above)
* Name of Reference (Outgroup) sequence (Must be the same as the reference sequence)
* Metadata file

``cladeomatic create --in_nwk examples/small_test/tree.nwk  --in_var examples/small_test/snps.vcf --in_meta examples/small_test/sample.meta.txt --outdir small_test_cladeomatic/ --root_name root.0 --reference examples/small_test/root.gbk``

Option 2 - Predefined groups

This mode will attempt to define a scheme based on a group manifest which meet membership size and SNP requirements.
Note every group id must be unique across all ranks to use this feature.Cladeomatic uses this for internal representation
of the genotypes but they can be mapped to whatever nomenclature is desired to have as an output

.. csv-table::
   :header: "sample_id", "invalid_genotype", "valid_genotype"
   :widths: 10,15,15

    "A", 0.1, 0.1
    "B", 0.1.1.1, 0.2.4.7
    "C", 0.1.1.2, 0.2.4.8
    "D", 0.1.1.2, 0.2.4.8
    "E", 0.2, 0.3
    "F", 0.2.1, 0.3.5
    "G", 0.2.2, 0.3.6
    "H", 0.2.2, 0.3.6

Input requirements are:

* TSV formatted group file (sample_id, genotype)
* VCF
* Reference sequence (.fasta / .gbk)
* Name of outgroup sequence
* Metadata file

``cladeomatic create --in_groups examples/small_test/groups.tsv --in_var examples/small_test/snps.vcf --in_meta examples/small_test/sample.meta.txt --outdir small_test_cladeomatic_groups/ --root_name root.0 --reference examples/small_test/root.gbk``

**Outputs:**

.. code-block:: console

    {Output folder name}
    ├── {prefix}-altseq.fasta - Artificial sequence which has a different base from the reference at every position in scheme
    ├── {prefix}-biohansel.fasta  - biohansel formatted kmer fasta file
    ├── {prefix}-biohansel.meta.txt - descriptions of biohansel kmers: kmername,target_position,target_base
    ├── {prefix}-clades.info.txt - Information on each individual clade, including supporting SNPs and metadata associations
    ├── {prefix}-dist.mat.txt - tab delimeted distance matrix from snp-dists
    ├── {prefix}-extracted.kmers.txt - Raw kmer output of extracted kmers with positions mapped
    ├── {prefix}-filtered.vcf - VCF file where invalid sites have been removed
    ├── {prefix}-genotypes.distance.txt - Histogram of node distances
    ├── {prefix}-genotypes.raw.txt - Tree or group file without filtering
    ├── {prefix}-genotypes.selected.txt - Nodes which meet the user criteria
    ├── {prefix}-genotypes.supported.txt - Nodes which were selected based on the supported nodes
    ├── {prefix}-kmer.scheme.txt - Cladeomatic kmer based scheme
    ├── {prefix}-params.log - Selected parameters for the run
    ├── {prefix}-sample.distances.html - Histogram of node distances
    ├── {prefix}-snps.scheme.txt - Cladeomatic SNP based scheme
    ├── {prefix}-snps.info.txt
    ├── pseudo.seqs.fasta - reconstructed fasta sequences based on reference sequence and vcf
    └──


**Genotype:**

Genotype samples using the developed scheme based on a VCF file with the same reference selected to build the scheme

Input requirements are:

* VCF
* Clade-O-Matic Scheme
* Metadata file (sample_id,genotype) * Produced by "create" {prefix}-genotypes.selected.txt
* (Optional) Metadata file (sample_id,genotype) * Produced by "create" {prefix}-genotypes.selected.txt

``cladeomatic genotype --in_var examples/small_test/snps.vcf --in_scheme examples/small_test/cladeomatic-snp.scheme.txt --sample_meta examples/small_test/sample.meta.txt --genotype_meta examples/small_test/genotype.meta.txt --outfile genotype.calls.txt``

VCF files will not include positions which are exclusively the reference sequence or missing and this poses an issue for calling
genotypes based on the VCF file where missing and reference state cannot be distinguished. A work around for this issue is the inclusion
of a sequence which is different from the reference sequence for every position targeted by the scheme. The create module generates a sequence
where every position used by the scheme is flipped to be a different base from the reference. This is not an ideal solution but it will allow
users to use the genotype module using SNIPPY-CORE with their query sequence and the "alt" sequence.


**Outputs:**

Outputs a file with the genotype calls for each input sample

**Benchmark Scheme:**

Benchmark the scheme based on the output of genotype tool. At this point only vcf based genotyping is supported

Input requirements are:

* TXT file produced by genotype module with predicted and expected genotypes, or tsv with predicted and submitted genotype information
* Clade-O-Matic scheme file used to call genotypes
* VCF
* Name of column for predicted genotype
* Name of column for submitted genotype

``cladeomatic benchmark --in_var examples/small_test/snps.vcf --in_scheme examples/small_test/cladeomatic-kmer.scheme.txt --in_genotype examples/small_test/genotype.calls.txt --submitted_genotype_col genotype --predicted_genotype_col predicted_genotype  --outdir benchmark``

The benchmark tool will identify the F1 scores for calling genotypes based on the provided scheme and will report per sample any sites which are responsible for
the submitted genotype not being called

**Outputs:**

.. code-block:: console

    OutputFolderName
    ├── {prefix}-scheme.scores.txt
    └── {prefix}-sample.results.txt


