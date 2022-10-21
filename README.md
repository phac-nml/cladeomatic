|logo|

|conda| |nbsp| |pypi| |nbsp|  |rtd| |nbsp| |license|


======  ===========
Master  |ci-master|
Dev     |ci-dev|
======  ===========

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


## Installation

Python dependencies (defined in the [requirements](https://github.com/phac-nml/cladeomatic/blob/main/requirements.txt) file, should be automatically installed when using conda or pip)

In addition to the python dependencies, Clade-o-Matic requires [Jellyfish 2.3.0](https://github.com/gmarcais/Jellyfish/)

Install the latest released version from conda:

        conda create -c bioconda -c conda-forge -n cladeomatic cladeomatic

Install using pip:

        pip install cladeomatic

Install the latest master branch version directly from Github:

        conda install jellyfish
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

    create  Define lineages and create a kmer scheme
    benchmark  Benchmark a kmer scheme
    test     Test parsityper functionality on a small dataset [ not implemented]
    version  Print version and exit

Quick start
=====
**Create scheme:**

Option 1 - De novo tree-based <br />
This mode will discover clades and lineages which meet membership size and SNP requirements. 
Input requirements are: 
* newick formatted tree
* VCF
* Reference sequence (.fasta / .gbk)
* Name of outgroup sequence
* Metadata file<br />
  

    cladeomatic create --in_nwk tree.nwk  --in_var variants.vcf --in_meta metadata.txt --outdir scheme/ --root_name ref --reference ref.gbk

Option 2 - Predefined groups <br />
This mode will attempt to define a scheme based on a group manifest which meet membership size and SNP requirements. 
Input requirements are: 
* TSV formatted group file (sample_id, genotype)
* VCF
* Reference sequence (.fasta / .gbk)
* Name of outgroup sequence
* Metadata file<br />
  

    cladeomatic create --in_groups groups.tsv --in_var variants.vcf --in_meta metadata.txt --outdir scheme/ --root_name ref --reference ref.gbk
  

**Outputs:**

```
OutputFolderName
├── {prefix}-clade.snp.histo.html [Tree Mode Only]
├── {prefix}-clades.info.txt
├── {prefix}-filt.kmers.txt
├── {prefix}-genotypes.raw.txt
├── {prefix}-genotypes.supported.txt
├── {prefix}-genotypes.selected.txt
├── {prefix}-genotype.consenus.fasta
├── {prefix}-jellyfish.counts.txt
├── {prefix}-scheme.txt
├── {prefix}-snps.all.txt
├── pseudo.seqs.fasta
├── samples.dists.matrix.csv [Tree Mode Only]
└──
```

**Benchmark Scheme:**
Benchmark the scheme using the original input VCF file and the set of genomes used to construct the scheme.
Input requirements are: 
* VCF
* Clade-O-Matic Scheme
* Metadata file (sample_id,genotype) * Produced by "create" {prefix}-genotypes.selected.txt
  

    cladeomatic benchmark --in_var variants.vcf --in_scheme cladeomatic-scheme.txt --in_meta metadata.txt --outdir benchmark/ 

Evaluate the results for any conflicting genotypes


**Outputs:**

```
OutputFolderName
├── {prefix}-scheme.scores.txt
└── {prefix}-scheme.calls.txt
```

## FAQ

## Citation

## Legal

Copyright Government of Canada 2022

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
