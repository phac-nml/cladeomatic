**Usage Guide**
===============

This section defines all the available options for running Clade-o-matic.

As a reminder, to get full help for a command use one of:

.. code-block:: console

    cladeomatic command -h
    cladeomatic command --help

**Main Options**

One of the main menu options below must be selected

.. csv-table::
   :header: "Argument", "Description"
   :widths: 15, 30

    "create", "Define lineages and create a kmer scheme"
    "gentoype", "Call genotypes from a VCF and a scheme file"
    "benchmark","Test developed scheme using labeled samples and scheme"
    "namer","Rename genotypes within a scheme"

**Create Scheme**

De novo tree data basic usage, ensure the ``--in_nwk`` option is selected:

``cladeomatic create --in_nwk examples/small_test/tree.nwk  --in_var examples/small_test/snps.vcf --in_meta examples/small_test/sample.meta.txt --outdir small_test_cladeomatic/ --root_name root.0 --reference examples/small_test/root.gbk``

Group data basic usage, ensure the ``--in_groups`` option is selected:

``cladeomatic create --in_groups examples/small_test/groups.tsv --in_var examples/small_test/snps.vcf --in_meta examples/small_test/sample.meta.txt --outdir small_test_cladeomatic_groups/ --root_name root.0 --reference examples/small_test/root.gbk``


.. csv-table:: **Create Scheme Options**
   :header: "Argument", "Required", "Description", "Input type", "Default Value"
   :widths: 30, 10, 40, 10, 10

    "\--in_var", True, "Either Variant Call SNP data (.vcf) or TSV SNP data (.txt)", String
    "\--in_nwk", False, "Newick Tree of strains", String
    "\--in_groups", False, "Tab delimited file of genotypes", String
    "\--in_meta", True, "Tab delimited file of metadata", String
    "\--reference", True, "Reference genbank or fasta sequence from VCF", String
    "\--outdir", True, "Output Directory to put results", String
    "\--prefix", False, "Prefix for output files", String, "cladeomatic"
    "\--root_name", False, "Name of sample to root tree", String
    "\--root_method", False, "Method to root tree (midpoint,outgroup)", String
    "\--klen", False, "kmer length", Integer, 18
    "\--min_members", False, "Minimum number of members for a clade to be valid", Integer, 1
    "\--min_snp_count", False, "Minimum number of unique SNPs for a clade to be valid", Integer, 1
    "\--max_snp_count", False, "Maximum number of SNPs to be selected for defining each genotype to prevent large numbers of redundant SNPs", Integer, 1
    "\--min_perc", False, "Minimum percentage of clade members to be positive for a kmer to be valid", Float, 0.1
    "\--max_site_ambig", False, "Maximum percentage of input sequences which can be missing a site for it to still be valid", Float, 0.25
    "\--max_states", False, "Maximum number of states for a position [A,T,C,G,N,-]", Integer, 6
    "\--max_ambig", False, "Maximum number of ambiguous bases allowed in a kmer", Integer, 0
    "\--rcor_thresh", False, "Correlation coefficient threshold", Float, 0.4
    "\--delim", False, "Genotype delimiter in group file", String
    "\--num_threads", False, "Number of threads to use", 1
    "\--no-plots", False, "Disable plotting",, True
    "\--keep_tmp", False, "Keep interim files",, True
    "\--debug", False, "Show debug information",, True
    "\--resume", False, "Resume previous analysis",, True
    "\--no_compression", False, "Skip compression of tree hierarchy",, True
    "\--force", False, "Force overwrite of existing results directory",, True
    "-V", False, "Provide the version number for this build"
    "\--version", False, "Provide the version number for this build"

**Genotyping**

Basic usage options:

``cladeomatic genotype --in_var examples/small_test/snps.vcf --in_scheme examples/small_test/cladeomatic-snp.scheme.txt --sample_meta examples/small_test/sample.meta.txt --genotype_meta examples/small_test/genotype.meta.txt --outfile genotype.calls.txt``

.. csv-table:: **Genotyping Options**
   :header: "Argument", "Required", "Description", "Input type", "Default Value"
   :widths: 30, 10, 40, 10, 10

    "\--in_var", True, "Either Variant Call SNP data (.vcf) or TSV SNP data (.txt)", String
    "\--in_scheme", True, "Tab delimited scheme file produced by clade-o-matic", String
    "\--sample_meta", True, "Tab delimited sample metadata", String
    "\--genotype_meta", False, "Tab delimited genotype metadata", String
    "\--outfile", True, "Output Directory to put results", String
    "\--max_missing_positions", False, "Maximum number of missing positions for the genotype", Integer, 1
    "\--num_threads", False, "Number of threads to use", Integer, 1
    "-V", False, "Provide the version number for this build"
    "\--version", False, "Provide the version number for this build"

**Benchmark Scheme:**

Basic usage options:

``cladeomatic benchmark --in_var examples/small_test/snps.vcf --in_scheme examples/small_test/cladeomatic-kmer.scheme.txt --in_genotype examples/small_test/genotype.calls.txt --submitted_genotype_col genotype --predicted_genotype_col predicted_genotype  --outdir benchmark``


.. csv-table:: **Benchmarking Options**
   :header: "Argument", "Required", "Description", "Input type", "Default Value"
   :widths: 30, 10, 40, 10, 10

    "\--in_genotype", True, "Genotype report made by genotyper", String
    "\--in_var", True, "Either Variant Call SNP data (.vcf) or TSV SNP data (.txt)", String
    "\--in_scheme", True, "Tab delimited scheme file produced by clade-o-matic", String
    "\--submitted_genotype_col", True, "Name of column containing submitted genotype", String
    "\--predicted_genotype_col", True, "Name of column containing predicted genotype", String
    "\--outdir", True, "Output Directory to put results", String
    "\--prefix", False, "Output Directory to put results", String, "cladeomatic"
    "\--debug", False, "Show debug information",, True
    "-V", False, "Provide the version number for this build"
    "\--version", False, "Provide the version number for this build"

**Namer**

Basic usage options:


.. csv-table:: **Namer Options**
   :header: "Argument", "Required", "Description", "Input type", "Default Value"
   :widths: 30, 10, 40, 10, 10

    "\--in_scheme", True, "Cladeomatic scheme file", String
    "\--in_name", True, "Tab delimited file of (node, name)", String
    "\--outfile", True, "Output file for updated scheme", String
    "-V", False, "Provide the version number for this build"
    "\--version", False, "Provide the version number for this build"