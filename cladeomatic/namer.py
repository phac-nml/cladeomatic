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


def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Clade-O-Matic: Genotyping scheme genotype namer v. {}".format(__version__),
        formatter_class=CustomFormatter)
    parser.add_argument('--in_scheme', type=str, required=True,
                        help='Cladeomatic scheme file')
    parser.add_argument('--in_names', type=str, required=True, help='Tab delimited file of genotype, name', default=None)
    parser.add_argument('--outfile', type=str, required=True, help='Output file for updated scheme')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()

def main():
    return