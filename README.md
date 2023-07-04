[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/cladeomatic/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/cladeomatic)
[![Conda](https://img.shields.io/conda/dn/bioconda/cladeomatic?color=green)](https://anaconda.org/bioconda/cladeomatic)
[![License: Apache-2.0](https://img.shields.io/github/license/phac-nml/cladeomatic)](https://www.apache.org/licenses/LICENSE-2.0)


<p align="left"><img src="logo.png" alt="Clade-o-matic" height="150" width="400"></p>

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Quick Start](#quick-start)
- [FAQ](#faq)
- [Citation](#citation)
- [Legal](#legal)
- [Contact](#contact)

## Introduction

Identification of population structure based on genomic data is a routine problem in molecular epidemiology and numerous approaches have been developed to partition large collections of sequences into clusters at different levels of granularity and resolution based on the goals of the researchers and underlying genetic diversity of the dataset. There is a confusing variety of terms in use to describe these units of genetic relatedness depending on the level of resolution provided and the specific taxon: clades, clones, clades, groups, genotypes, haplotypes, lineages, and subtypes.

Clade-o-matic is a phylogenetic approach for identification of hierarchal genotypes based on canonical SNPs that are exclusive and conserved within phylogenetic clades. It works best with organisms that have low rates of recombination and limited genetic diversity, but it is possible to analyze species-wide phylogenies. A major distinction of Clade-o-matic from other population structure tools such as rheirbaps is that it provides the distinguishing features for each clade in the form of a SNP typing scheme with unique kmers for identification of these SNPs from raw sequencing data. Furthermore, Clade-o-matic does not use a fixed set of levels for hierarchy as there will be uneven groupings between different clades. It uses the underlying distance distributions to compress the hierarchy into the minimum set of levels that retain the tree structure. The granularity of the scheme can be changed by specifiying different numbers of SNPs required to support a clade along with the number of members assigned to a clade for it to be considered valid.

### Check the [documentation](https://cladeomatic-documentation.readthedocs.io/en/latest/) for implementation details and guidance on using Clade-O-Matic.

## Installation

Python dependencies (defined in the [requirements](https://github.com/phac-nml/cladeomatic/blob/main/requirements.txt) file, should be automatically installed when using conda or pip)

In addition to the python dependencies, Clade-o-Matic requires:
[Jellyfish 2.3.0](https://github.com/gmarcais/Jellyfish/)
[snp-dists 0.8.2](https://github.com/tseemann/snp-dists/)


Install the latest released version from conda:

        conda create -c bioconda -c conda-forge -n cladeomatic cladeomatic

Install using pip:

        pip install cladeomatic

Install the latest master branch version directly from Github:

        conda install jellyfish snp-dists
        pip install git+https://github.com/phac-nml/cladeomatic.git



## Usage
If you run ``cladeomatic``, you should see the following usage statement:

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
    
Quick start
=====
**Create scheme:**

Option 1 - De novo tree-based <br />
This mode will discover clades and lineages which meet membership size and SNP requirements. 
Input requirements are: 
* newick formatted tree
* Reference (Outgroup) sequence (.fasta / .gbk)
* VCF (Must use the same reference sequence as above)
* Name of Reference (Outgroup) sequence (Must be the same as the reference sequence)
* Metadata file<br />


        cladeomatic create --in_nwk examples/small_test/tree.nwk  --in_var examples/small_test/snps.vcf --in_meta examples/small_test/sample.meta.txt --outdir small_test_cladeomatic/ --root_name root.0 --reference examples/small_test/root.gbk

Option 2 - Predefined groups <br />
This mode will attempt to define a scheme based on a group manifest which meet membership size and SNP requirements.
Note every group id must be unique accross all ranks to use this feature.Cladeomatic uses this for internal representation 
of the genotypes but they can be mapped to whatever nomenclature is desired to have as an output

| sample_id     | invalid_genotype | valid_genotype |
| ----------- | ----------- |----------- | 
| A      | 0.1      |   0.1         |
| B   | 0.1.1.1        |     0.2.4.7       |
| C     | 0.1.1.2       |    0.2.4.8        |
| D   | 0.1.1.2        |       0.2.4.8       |
| E      | 0.2      |       0.3     |
| F   | 0.2.1        |      0.3.5     |
| G   | 0.2.2        |       0.3.6   |
| H   | 0.2.2        |      0.3.6      |

Input requirements are: 
* TSV formatted group file (sample_id, genotype)
* VCF
* Reference sequence (.fasta / .gbk)
* Name of outgroup sequence
* Metadata file<br />


        cladeomatic create --in_groups examples/small_test/groups.tsv --in_var examples/small_test/snps.vcf --in_meta examples/small_test/sample.meta.txt --outdir small_test_cladeomatic_groups/ --root_name root.0 --reference examples/small_test/root.gbk
  

**Outputs:**

```
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
```

**Genotype:**

Genotype samples using the developed scheme based on a VCF file with the same reference selected to build the scheme
Input requirements are: 
* VCF
* Clade-O-Matic Scheme
* Metadata file (sample_id,genotype) * Produced by "create" {prefix}-genotypes.selected.txt
* (Optional) Metadata file (sample_id,genotype) * Produced by "create" {prefix}-genotypes.selected.txt


    cladeomatic genotype --in_var examples/small_test/snps.vcf --in_scheme examples/small_test/cladeomatic-snp.scheme.txt --sample_meta examples/small_test/sample.meta.txt --genotype_meta examples/small_test/genotype.meta.txt --outfile genotype.calls.txt

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


    cladeomatic benchmark --in_var examples/small_test/snps.vcf --in_scheme examples/small_test/cladeomatic-kmer.scheme.txt --in_genotype examples/small_test/genotype.calls.txt --submitted_genotype_col genotype --predicted_genotype_col predicted_genotype  --outdir benchmark

The benchmark tool will identify the F1 scores for calling genotypes based on the provided scheme and will report per sample any sites which are responsible for 
the submitted genotype not being called

**Outputs:**

```
OutputFolderName
├── {prefix}-scheme.scores.txt
└── {prefix}-sample.results.txt
```

## FAQ

## Citation

## Legal

Copyright Government of Canada 2023

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.


## Contact

**James Robertson**: james.robertson@phac-aspc.gc.ca
