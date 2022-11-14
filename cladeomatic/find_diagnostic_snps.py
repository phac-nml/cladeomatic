import copy
import statistics
import time
import random
import logging
import os
import copy
import psutil
import shutil
import math
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from statistics import mean
import numpy as np
import pandas as pd
import ray
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.signal import find_peaks
from scipy.stats import spearmanr, pearsonr

from cladeomatic.constants import SCHEME_HEADER, IUPAC_LOOK_UP, CLADE_DEFINING_SNPS_HEADER
from cladeomatic.utils import init_console_logger, calc_ARI, calc_AMI, calc_shanon_entropy, fisher_exact
from cladeomatic.utils import parse_metadata
from cladeomatic.utils.jellyfish import run_jellyfish_count, parse_jellyfish_counts
from cladeomatic.utils.kmerSearch import SeqSearchController, revcomp
from cladeomatic.utils.phylo_tree import parse_tree, prune_tree, tree_to_distance_matrix, \
    get_pairwise_distances_from_matrix
from cladeomatic.utils.seqdata import create_pseudoseqs_from_vcf, parse_reference_gbk, \
    create_aln_pos_from_unalign_pos_lookup, calc_homopolymers
from cladeomatic.utils.vcfhelper import vcfReader
from cladeomatic.utils.visualization import plot_bar
from cladeomatic.utils.seqdata import get_variants
from cladeomatic.writers import write_snp_report, write_kmers, write_genotypes, write_clade_snp_report, write_node_report, print_params
from cladeomatic.version import __version__

fasta_file = '/Users/jrobertson/Desktop/__2022-Cladeomatic/aroC/aroC.aligned.fas'
genotypes_file = '/Users/jrobertson/Desktop/__2022-Cladeomatic/aroC/cladeomatic/cladeomatic-genotypes.raw.txt'
genotypes_df = pd.read_csv(genotypes_file,header=0,sep="\t")
genotype_map = {}
genotypes = set()
for row in genotypes_df.itertuples():
    sample_id = row.sample_id
    genotype = row.genotype
    genotype_map[sample_id] = genotype
    genotypes.add(genotype)


consensus = []
genotype_counts = {}
bases = set(['A','T','C','G','-','N'])
with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        id = str(record.id)
        seq = str(record.seq).upper()
        seq_len = len(seq)
        genotype = genotype_map[sample_id]
        if len(consensus) == 0:
            for i in range(0,seq_len):
                consensus.append({'A':0,'T':0,'C':0,'G':0,'N':0,'-':0})
        if not genotype in genotype_counts:
            genotype_counts[genotype] = []
            for i in range(0,seq_len):
                genotype_counts[genotype].append({'A':0,'T':0,'C':0,'G':0,'N':0,'-':0})
            nodes = genotype.split('.')
            for n in nodes:
                genotype_counts[n] = []
                for i in range(0, seq_len):
                    genotype_counts[n].append({'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, '-': 0})

        for i in range(0,seq_len):
            base = seq[i]
            if base not in bases:
                base = 'N'
            consensus[i][base]+=1
            genotype_counts[genotype][i][base]+=1
            for n in nodes:
                genotype_counts[n][i][base]+=1
handle.close()

valid_bases = ['A','T','C','G']
for genotype in genotype_counts:
    for i in range(0,seq_len):
        total = 0
        for b in valid_bases:
            total +=genotype_counts[genotype][i][b]
        nt = max(genotype_counts[genotype][i], key=genotype_counts[genotype][i].get)
        nt_count = genotype_counts[genotype][i][nt]
        #print("{}\t{}\t{}\t{}\t{}\t{}".format(genotype, i, nt,nt_count,total,consensus[i][nt]))
        if nt == 'N' or nt == '-' or total != nt_count:
            continue
        consensus_count = consensus[i][nt]
        if consensus_count != nt_count:
            continue
        print("{}\t{}\t{}".format(genotype,i,nt))






