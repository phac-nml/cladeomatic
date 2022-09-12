import copy
import tempfile
import shutil
import logging
import math
import os
import psutil
import shutil
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

from cladeomatic.constants import SCHEME_HEADER
from cladeomatic.utils import init_console_logger, calc_ARI, calc_AMI, calc_shanon_entropy, fisher_exact
from cladeomatic.utils import parse_metadata
from cladeomatic.utils.jellyfish import run_jellyfish_count, parse_jellyfish_counts
from cladeomatic.utils.kmerSearch import SeqSearchController
from cladeomatic.utils.phylo_tree import parse_tree, prune_tree, tree_to_distance_matrix, \
    get_pairwise_distances_from_matrix, annotate_tree, plot_single_rep_tree
from cladeomatic.utils.seqdata import create_pseudoseqs_from_vcf, parse_reference_gbk, \
    create_aln_pos_from_unalign_pos_lookup, calc_homopolymers
from cladeomatic.utils.vcfhelper import vcfReader
from cladeomatic.utils.visualization import plot_bar
from cladeomatic.utils.seqdata import get_variants
from cladeomatic.version import __version__


def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Clade-O-Matic: Genotyping scheme development v. {}".format(__version__),
        formatter_class=CustomFormatter)
    parser.add_argument('--in_var', type=str, required=True,
                        help='Either Variant Call SNP data (.vcf) or TSV SNP data (.txt)')
    parser.add_argument('--in_nwk', type=str, required=True, help='Newick Tree of strains')
    parser.add_argument('--in_meta', type=str, required=True, help='Tab delimited file of metadata', default=None)
    parser.add_argument('--reference', type=str, required=True, help='Reference genbank or fasta sequence from VCF')
    parser.add_argument('--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('--prefix', type=str, required=False, help='Prefix for output files', default='cladeomatic')
    parser.add_argument('--root_name', type=str, required=False, help='Name of sample to root tree', default='')
    parser.add_argument('--root_method', type=str, required=False, help='Method to root tree (midpoint,outgroup)',
                        default=None)
    parser.add_argument('--klen', type=int, required=False, help='kmer length', default=18)
    parser.add_argument('--min_members', type=int, required=False,
                        help='Minimum number of members for a clade to be valid', default=1)
    parser.add_argument('--min_snp_count', type=int, required=False,
                        help='Minimum number of unique snps for a clade to be valid',
                        default=1)
    parser.add_argument('--min_perc', type=float, required=False,
                        help='Minimum percentage of clade members to be positive for a kmer to be valid', default=1)
    parser.add_argument('--max_states', type=int, required=False,
                        help='Maximum number of states for a position [A,T,C,G,N,-]',
                        default=6)
    parser.add_argument('--rcor_thresh', type=float, required=False, help='Correlation coefficient threshold',
                        default=0.4)
    parser.add_argument('--num_threads', type=int, required=False, help='Number of threads to use', default=1)
    parser.add_argument('--no_plots', required=False, help='Disable plotting', action='store_true')
    parser.add_argument('--keep_tmp', required=False, help='Keep interim files', action='store_true')
    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()


def generate_genotypes(tree):
    '''
    Accepts a ETE3 tree object and returns the heirarchal node id's for every leaf
    :param tree: ETE3 tree object
    :return: dict of leaf names and list of node heirarchy
    '''
    samples = tree.get_leaf_names()
    geno = {}
    for sample_id in samples:
        geno[sample_id] = []
    for node in tree.traverse("preorder"):
        node_id = node.name
        if node.is_leaf():
            continue
        leaves = node.get_leaf_names()
        for sample_id in leaves:
            geno[sample_id].append(node_id)
    return geno


def write_genotypes(genotypes, outfile, delimeter='.'):
    '''
    Accepts a list of sample genotype heirarchies and writes them to a file
    :param genotypes: dict of leaf names with heirarchal list of node id's to represent tree structure
    :param outfile: str
    :param delimeter: str
    :return:
    '''
    fh = open(outfile, 'w')
    fh.write("sample_id\tgenotype\n")
    for sample_id in genotypes:
        genotype = delimeter.join(genotypes[sample_id])
        fh.write("{}\t{}\n".format(sample_id, genotype))
    fh.close()
    return


def write_snp_report(snp_data, outfile):
    '''
    Accepts snp_data dict data structure and writes snp details
    :param snp_data: dict() {[chrom_id] : {[position]: {[base]: :dict()}}}
    :param outfile: str
    :return:
    '''
    fh = open(outfile, 'w')
    fh.write("chrom\tpos\tbase\tclade_id\tis_canonical\tnum_clade_members\tnum_members\tis_ref\n")
    for chrom in snp_data:
        for pos in snp_data[chrom]:
            for base in snp_data[chrom][pos]:
                row = [chrom, pos, base,
                       snp_data[chrom][pos][base]['clade_id'],
                       snp_data[chrom][pos][base]['is_canonical'],
                       snp_data[chrom][pos][base]['num_clade_members'],
                       snp_data[chrom][pos][base]['num_members'],
                       snp_data[chrom][pos][base]['is_ref']]
                fh.write("{}\n".format("\t".join([str(x) for x in row])))
    fh.close()
    return


def write_node_report(clade_data, outfile):
    '''
    Writes the Node information data to a tsv file
    :param clade_data: dict Node information structure
    :param outfile: str
    :return:
    '''
    fh = open(outfile, 'w')
    fh.write("clade_id\tchrom\tpos\tbase\tis_canonical\tnum_clade_members\tnum_members\tis_ref\n")
    for clade_id in clade_data:
        row = [clade_id]
        for feat in clade_data[clade_id]:
            value = clade_data[clade_id][feat]
            if isinstance(value, list):
                value = ";".join([str(x) for x in value])
            row.append(value)
        fh.write("{}\n".format("\t".join([str(x) for x in row])))
    fh.close()

    return


def identify_canonical_snps(ete_tree_obj, vcf_file, min_member_count=1):
    '''
    Accepts SNP data and tree to identify which SNPs correspond to a specific node on the tree
    :param ete_tree_obj: ETE3 tree object
    :param vcf_file: str path to vcf or tsv snp data
    :return: dict of snp_data data structure
    '''
    vcf = vcfReader(vcf_file)
    data = vcf.process_row()
    samples = vcf.samples
    num_samples = len(samples)
    snps = {}
    count_snps = 0
    if data is not None:
        count_snps += 1
    node_cache = ete_tree_obj.get_cached_content()
    node_lookup = {}
    for node_id in ete_tree_obj.traverse():
        node_members = list(node_cache[node_id])
        node_lookup[node_id.name] = set()
        for leaf in node_members:
            node_lookup[node_id.name].add(leaf.name)

    while data is not None:
        chrom = data['#CHROM']
        pos = int(data['POS']) - 1
        if not chrom in snps:
            snps[chrom] = {}
        snps[chrom][pos] = {}
        assignments = {}
        for sample_id in samples:
            base = data[sample_id]
            if not base in assignments:
                assignments[base] = []
            assignments[base].append(sample_id)
        is_ambig = False
        if 'N' in assignments:
            is_ambig = True
        for base in assignments:
            if not base in ['A', 'T', 'C', 'G']:
                continue
            in_samples = set(assignments[base])
            is_ref = base == data['REF']
            is_canonical = False
            best_node_id = ''
            best_count = num_samples

            for clade_id in node_lookup:
                node_leaves = node_lookup[clade_id]
                if is_ambig:
                    node_leaves = set(node_leaves) - set(assignments['N'])
                i1 = (in_samples & node_leaves)
                i2 = (node_leaves & in_samples)
                num_out_members = len(node_leaves - i2)
                if len(in_samples - node_leaves) == 0 & num_out_members < best_count:
                    best_count = num_out_members
                    best_node_id = clade_id
                if num_out_members == 0 and len(in_samples) == len(i1):
                    is_canonical = True
                    break
            clade_id = best_node_id
            clade_total_members = len(node_lookup[clade_id])

            snps[chrom][pos][base] = {'clade_id': clade_id, 'is_canonical': is_canonical,
                                      'num_clade_members': clade_total_members, 'num_members': clade_total_members,
                                      'is_ref': is_ref}
        count_snps += 1
        data = vcf.process_row()

    return snps


def snp_based_filter(snp_data, min_count, max_states):
    '''
    Accepts a snp_data dict and filters the dictionary to return only those snps which are present in >=min_count samples
    and have no more than max_state bases at a given position
    :param snp_data: dict() {[chrom_id] : {[position]: {[base]: :dict()}}}
    :param min_count: int
    :param max_states: int
    :return: dict
    '''
    filtered = {}
    for chrom in snp_data:
        filtered[chrom] = {}
        for pos in snp_data[chrom]:
            for base in snp_data[chrom][pos]:
                if snp_data[chrom][pos][base]['is_canonical']:
                    if snp_data[chrom][pos][base]['num_members'] >= min_count:
                        if not pos in filtered[chrom]:
                            filtered[chrom][pos] = {}
                        filtered[chrom][pos][base] = snp_data[chrom][pos][base]
            if pos in filtered[chrom] and len(filtered[chrom][pos]) > max_states:
                del (filtered[chrom][pos])

    return filtered


def get_valid_nodes(snp_data, min_snp_count):
    '''
    Accepts snp_data structure and identifies which nodes have sufficient canonical SNPS to consider the node
    valid
    :param snp_data: dict() {[chrom_id] : {[position]: {[base]: :dict()}}}
    :param min_snp_count: int
    :return: clade data structure of valid nodes in the tree
    '''
    nodes = {}
    for chrom in snp_data:
        for pos in snp_data[chrom]:
            for base in snp_data[chrom][pos]:
                if snp_data[chrom][pos][base]['is_canonical']:
                    clade_id = snp_data[chrom][pos][base]['clade_id']

                    if not clade_id in nodes:
                        nodes[clade_id] = {
                            'chr': [],
                            'pos': [],
                            'bases': [],
                            'min_dist': 0,
                            'max_dist': 0,
                            'ave_dist': 0,
                            'entropies': {},
                            'fisher': {},
                            'ari': {},
                            'ami': {}
                        }
                    nodes[clade_id]['chr'].append(chrom)
                    nodes[clade_id]['pos'].append(pos)
                    nodes[clade_id]['bases'].append(base)
    filtered = {'0': {'chr': [],
                      'pos': [],
                      'bases': [],
                      'min_dist': 0,
                      'max_dist': 0,
                      'ave_dist': 0,
                      'entropies': {},
                      'fisher': {},
                      'ari': {},
                      'ami': {}}}

    for clade_id in nodes:
        if len(nodes[clade_id]['pos']) >= min_snp_count:
            filtered[clade_id] = nodes[clade_id]

    return filtered


def calc_node_distances(clade_data, ete_tree_obj):
    '''
    Accepts a dict of clade_data and for each node calculates the distance of all leaves to the node
    :param clade_data: dict clade data structure of tree nodes
    :param ete_tree_obj: ETE3 tree obj
    :return: dict filtered clade data structure
    '''
    node_cache = ete_tree_obj.get_cached_content()

    for node_id in ete_tree_obj.traverse():
        clade_id = node_id.name
        if clade_id not in clade_data:
            continue
        node_members = list(node_cache[node_id])
        distances = []
        for leaf in node_members:
            dist = node_id.get_distance(leaf)
            distances.append(dist)
        min_dist = 0
        max_dist = 0
        ave_dist = 0
        if len(distances) > 0:
            min_dist = min(distances)
            max_dist = max(distances)
            ave_dist = mean(distances)
        clade_data[clade_id]['min_dist'] = min_dist
        clade_data[clade_id]['max_dist'] = max_dist
        clade_data[clade_id]['ave_dist'] = ave_dist

    return clade_data


def calc_node_associations(metadata, clade_data, ete_tree_obj):
    '''
    :param metadata:
    :param clade_data:
    :param ete_tree_obj:
    :return:
    '''
    samples = set(metadata.keys())
    num_samples = len(samples)
    for clade_id in clade_data:
        node = ete_tree_obj.search_nodes(name=clade_id)[0]
        in_members = set(node.get_leaf_names())
        features = {}
        genotype_assignments = []
        metadata_labels = {}
        ftest = {}
        metadata_counts = {}
        for sample_id in samples:
            if sample_id in in_members:
                genotype_assignments.append(1)
            else:
                genotype_assignments.append(0)
            for field_id in metadata[sample_id]:
                value = metadata[sample_id][field_id]
                if not field_id in metadata_labels:
                    metadata_counts[field_id] = {}
                    metadata_labels[field_id] = []
                    ftest[field_id] = {}
                if not value in ftest[field_id]:
                    metadata_counts[field_id][value] = 0
                    ftest[field_id][value] = {
                        'pos-pos': set(),
                        'pos-neg': set(),
                        'neg-pos': set(),
                        'neg-neg': set()
                    }
                metadata_counts[field_id][value] += 1
                if sample_id in in_members:
                    ftest[field_id][value]['pos-pos'].add(sample_id)
                else:
                    ftest[field_id][value]['neg-pos'].add(sample_id)

                metadata_labels[field_id].append(value)
                if not field_id in features:
                    features[field_id] = {}

                if not value in features[field_id]:
                    features[field_id][value] = 0
                if sample_id in in_members:
                    features[field_id][value] += 1

        for field_id in metadata_labels:
            category_1 = []
            category_2 = []
            for idx, value in enumerate(metadata_labels[field_id]):
                if isinstance(value,float):
                    if math.isnan(value):
                        continue
                value = str(value)
                if len(value) == 0 or value == 'nan':
                    continue
                category_1.append(genotype_assignments[idx])
                category_2.append(metadata_counts[field_id][value])

            clade_data[clade_id]['ari'][field_id] = calc_ARI(category_1, category_2)
            clade_data[clade_id]['ami'][field_id] = calc_AMI(category_1, category_2)
            clade_data[clade_id]['entropies'][field_id] = calc_shanon_entropy(category_2)
            clade_data[clade_id]['fisher'][field_id] = {}

            for value in ftest[field_id]:
                ftest[field_id][value]['neg-neg'] = (
                            ftest[field_id][value]['pos-pos'] | ftest[field_id][value]['neg-pos'])
                ftest[field_id][value]['pos-neg'] = in_members - ftest[field_id][value]['neg-neg']
                table = [
                    [len(ftest[field_id][value]['pos-pos']),
                     len(ftest[field_id][value]['neg-neg'] | ftest[field_id][value]['pos-neg'])],
                    [len(ftest[field_id][value]['neg-pos'] | ftest[field_id][value]['pos-pos']),
                     len(ftest[field_id][value]['neg-neg'] | ftest[field_id][value]['neg-pos'])]
                ]
                oddsr, p = fisher_exact(table, alternative='greater')
                clade_data[clade_id]['fisher'][field_id][value] = {'oddsr': oddsr, 'p': p}

    return clade_data





def sample_dist_summary(distances, num_bins=25):
    '''
    Divides the distances into equal sized bins
    :param distances: list of distances between nodes
    :return: dict
    '''
    df = pd.DataFrame(distances, columns=['dist'])
    qcut_series, qcut_intervals = pd.qcut(df['dist'], q=num_bins, retbins=True, duplicates='drop')
    tmp = pd.DataFrame(qcut_series.value_counts(sort=True))
    tmp = tmp.sort_index()
    tmp.reset_index(inplace=True)
    tmp = tmp.rename(columns={"index": "interval", "dist": "count"})
    return {'summary': df.describe().to_dict(), 'results': tmp, 'intervals': qcut_intervals}


def find_dist_peaks(values):
    '''
    :param values: list
    :return:
    '''
    return (find_peaks(np.array(values)))


def find_dist_troughs(values):
    '''
    :param values: list
    :return:
    '''
    return (find_peaks(-np.array(values)))


def node_list_attr(clade_data, attribute):
    '''
    Parameters
    ----------
    clade_data
    attribute
    Returns
    -------
    '''
    values = []
    for clade_id in clade_data:
        value = clade_data[clade_id][attribute]
        if isinstance(value, list):
            value = len(value)
        values.append(value)
    return values


def temporal_signal(metadata, clade_data, ete_tree_obj, rcor_thresh):
    '''
    Parameters
    ----------
    metadata
    clade_data
    ete_tree_obj
    rcor_thresh
    Returns
    -------
    '''
    for clade_id in clade_data:
        node = ete_tree_obj.search_nodes(name=clade_id)[0]
        in_members = set(node.get_leaf_names())
        distances = []
        years = []
        for sample_id in in_members:
            value = metadata[sample_id]['year']
            if value == 'nan' or len(value) == 0:
                continue
            dist = node.get_distance(sample_id)
            distances.append(dist)

            years.append(int(float(metadata[sample_id]['year'])))

        R = 0
        P = 1
        Rpear = 0
        Ppear = 1
        if len(np.asarray(years)) >= 2 and len(np.asarray(distances)) >= 2:
            R, P = spearmanr(np.asarray(years), np.asarray(distances))
            Rpear, Ppear = pearsonr(np.asarray(years), np.asarray(distances))
        clade_data[clade_id]['spearmanr'] = R
        clade_data[clade_id]['spearmanr_pvalue'] = P
        clade_data[clade_id]['pearsonr'] = Rpear
        clade_data[clade_id]['pearsonr_pvalue'] = Ppear
        is_present = False
        if clade_data[clade_id]['spearmanr'] > rcor_thresh or clade_data[clade_id]['pearsonr'] > rcor_thresh:
            is_present = True
        clade_data[clade_id]['is_temporal_signal_present'] = is_present
    return clade_data


def as_range(values):
    '''
    Parameters
    ----------
    values list of integers
    Returns list of tuples with the start and end index of each range
    -------
    '''
    values = sorted(values)
    num_val = len(values)
    if num_val == 0:
        return []
    if num_val == 1:
        return [(0, 0)]
    ranges = []
    i = 0
    while i < num_val - 1:

        s = i
        for k in range(i + 1, num_val):
            if values[k] - values[i] != 1:
                i += 1
                break
            i += 1
        if values[k] - values[i] > 1:
            ranges.append((s, k - 1))
        else:
            if k == num_val - 1:
                ranges.append((s, k))
            else:
                ranges.append((s, k - 1))
    return ranges


def get_nomenclature_ranges(counts):
    '''
    Parameters
    ----------
    counts: list of frequencies
    Returns list of tuples of compressed ranges
    -------
    '''
    (troughs, _) = find_dist_troughs(counts)
    troughs = list(troughs)
    min_value = min(counts)
    for i, value in enumerate(counts):
        if value == min_value:
            troughs.append(i)
    troughs = sorted(list(set(troughs)))
    c = as_range(troughs)
    ranges = []
    for idx, tuple in enumerate(c):
        ranges.append((troughs[tuple[0]], troughs[tuple[1]]))

    return ranges


def select_nodes(clade_data, dist_ranges, min_count=1):
    '''
    Parameters
    ----------
    clade_data
    dist_ranges
    min_count
    Returns
    -------
    '''
    num_ranks = len(dist_ranges)
    candidate_nodes = []
    for i in range(0, num_ranks):
        candidate_nodes.append([])
    candidate_nodes[-1].append('0')
    for clade_id in clade_data:

        if len(clade_data[clade_id]['chr']) < min_count:
            continue
        ave_dist = clade_data[clade_id]['ave_dist']

        for i, tuple in enumerate(dist_ranges):
            s = tuple[0]
            e = tuple[1]
            if ave_dist >= s and ave_dist < e:
                candidate_nodes[i].append(clade_id)

    return candidate_nodes


def select_bifucating_nodes(ete_tree_obj, clade_data):
    '''
    Parameters
    ----------
    ete_tree_obj
    clade_data
    Returns
    -------
    '''
    bifurcating_nodes = []
    for clade_id in clade_data:
        children = ete_tree_obj.search_nodes(name=clade_id)[0].get_children()
        node_count = 0
        if clade_id in bifurcating_nodes:
            continue
        nodes = [clade_id]
        for child in children:
            if not child.is_leaf() and child.name in clade_data:
                node_count += 1
                nodes.append(child.name)
        if node_count > 1:
            bifurcating_nodes += nodes
    return bifurcating_nodes


def create_compressed_hierarchy(ete_tree_obj, selected_nodes):
    '''
    Parameters
    ----------
    ete_tree_obj
    selected_nodes
    Returns
    -------
    '''
    selected_nodes.reverse()
    sample_genotypes = generate_genotypes(ete_tree_obj)
    valid_nodes = ['0']
    filt_nodes = []
    for idx, nodes in enumerate(selected_nodes):
        for ni, node_id in enumerate(nodes):
            valid_nodes.append(node_id)
        if len(nodes) > 0:
            filt_nodes.append(nodes)


    for sample_id in sample_genotypes:

        genotype = sample_genotypes[sample_id]
        filt = []
        for node_id in genotype:
            if node_id in valid_nodes:
                filt.append(node_id)

        for idx, nodes in enumerate(filt_nodes):
            ovl = set(nodes) & set(filt)
            if len(ovl) > 1:
                match = False
                genotype = []
                for i in range(0, len(filt)):
                    node = filt[i]
                    if node in ovl:
                        if not match:
                            genotype.append(node)
                            match = True
                            continue
                    else:
                        genotype.append(node)
                filt = genotype
        sample_genotypes[sample_id] = filt
    valid_nodes = set(['0'])
    terminal_nodes = set()
    for sample_id in sample_genotypes:
        valid_nodes = valid_nodes & set(sample_genotypes[sample_id])
        terminal_nodes.add(sample_genotypes[sample_id][-1])
    valid_nodes = valid_nodes | terminal_nodes
    return valid_nodes


def validate_file(file):
    '''
    Parameters
    ----------
    file
    Returns
    -------
    '''
    if os.path.isfile(file) and os.path.getsize(file) > 9:
        return True
    else:
        return False


def isint(str):
    '''
    Parameters
    ----------
    str
    Returns
    -------
    '''
    try:
        int(str)
        return True
    except ValueError:
        return False


def add_kmer_positions(kmer_mapping_info, klen, pseudo_seq_file):
    '''
    Parameters
    ----------
    kmer_mapping_info
    klen
    pseudo_seq_file
    Returns
    -------
    '''
    seqs_to_check = {}

    for kIndex in kmer_mapping_info:
        match_index = int(kmer_mapping_info[kIndex]['match_index'])
        seq_id = kmer_mapping_info[kIndex]['seq_id']
        kmer_mapping_info[kIndex]['uStart'] = match_index - klen + 1
        kmer_mapping_info[kIndex]['uEnd'] = match_index
        if not seq_id in seqs_to_check:
            seqs_to_check[seq_id] = []
        seqs_to_check[seq_id].append(kIndex)
    with open(pseudo_seq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = str(record.id)
            if id not in seqs_to_check:
                continue
            seq = str(record.seq).upper().replace('-', '')
            aln_lookup = create_aln_pos_from_unalign_pos_lookup(seq)
            for kIndex in seqs_to_check[id]:
                uStart = kmer_mapping_info[kIndex]['uStart']
                uEnd = kmer_mapping_info[kIndex]['uEnd']
                kmer_mapping_info[kIndex]['aStart'] = aln_lookup[uStart]
                kmer_mapping_info[kIndex]['aEnd'] = aln_lookup[uEnd]
        handle.close()
    return kmer_mapping_info


def group_kmers_by_pos(kmer_mapping_info):
    '''
    Parameters
    ----------
    kmer_mapping_info
    Returns
    -------
    '''
    groups = {}
    id = 0
    for kIndex in kmer_mapping_info:
        aStart = kmer_mapping_info[kIndex]['aStart']
        if not aStart in groups:
            groups[aStart] = {}
        groups[aStart][id] = kmer_mapping_info[kIndex]
        id += 1
    return groups


def select_kmers(kmer_groups, clade_info, klen, min_kmers=1, max_kmers=100):
    '''
    Parameters
    ----------
    kmer_groups
    clade_info
    klen
    min_kmers
    max_kmers
    Returns
    -------
    '''
    selected_kmers = {}
    variant_positions = set()
    conserved_positions = set()
    for pos in kmer_groups:
        num_kmers = kmer_groups[pos]['num_kmers']
        if num_kmers == 1:
            conserved_positions.add(pos)
    conserved_positions = sorted(list(conserved_positions))
    conserved_range_indexes = as_range(conserved_positions)
    conserved_ranges = []
    for i in range(0, len(conserved_range_indexes)):
        s = conserved_range_indexes[i][0]
        e = conserved_range_indexes[i][1]
        conserved_ranges.append((conserved_positions[s], conserved_positions[e]))

    positions = clade_info['0']['pos']
    for i in range(0,len(conserved_ranges)):
        pos = conserved_ranges[i][0]
        if len(kmer_groups[pos]['kmers']) > 1:
            continue
        positions.append(pos)
        clade_info['0']['chr'].append('1')
        kIndex = list(kmer_groups[pos]['kmers'].keys())[0]

        kseq = kmer_groups[pos]['kmers'][kIndex]
        base = kseq[0]
        clade_info['0']['bases'].append(base)
        selected_kmers[pos] = {
            'out_group': {},
            'in_group': {base: [kIndex]},
        }

    clade_info['0']['pos'] = positions

    for clade_id in clade_info:
        if clade_id == '0':
            continue
        positions = clade_info[clade_id]['pos']
        variant_positions = set(positions) | variant_positions
        bases = clade_info[clade_id]['bases']
        in_group = []
        out_group = []
        for idx, pos in enumerate(positions):
            base = bases[idx]
            s = pos - klen + 1
            if s < 0:
                s = 0
            e = pos
            group_ids = range(s, e)
            for group_id in group_ids:
                if group_id not in kmer_groups:
                    continue
                in_group = {}
                out_group = {}
                for kIndex in kmer_groups[group_id]['kmers']:
                    kseq = kmer_groups[group_id]['kmers'][kIndex]
                    bpos = pos - group_id
                    if kseq[bpos] == base:
                        if not base in in_group:
                            in_group[kseq[bpos]] = []
                        in_group[kseq[bpos]].append(kIndex)
                    else:
                        if not kseq[bpos] in out_group:
                            out_group[kseq[bpos]] = []

                        out_group[kseq[bpos]].append(kIndex)
                if len(in_group) > 0 and len(out_group) > 0:
                    break
            if len(in_group) > 0 and len(out_group) > 0:
                selected_kmers[pos] = {
                    'out_group': out_group,
                    'in_group': in_group,
                }


    group_map = {}
    for group_id in kmer_groups:
        for kIndex in kmer_groups[group_id]['kmers']:
            group_map[kIndex] = group_id

    selected_kmers = minimize_kmers(selected_kmers, kmer_groups, clade_info, klen, min_kmers, max_kmers)

    kmers = {}
    for pos in selected_kmers:
        kmers[pos] = {}
        for group in selected_kmers[pos]:
            for base in selected_kmers[pos][group]:
                kmer_ids = selected_kmers[pos][group][base]
                kmers[pos][base] = {}
                for kIndex in kmer_ids:
                    group_id = group_map[int(kIndex)]
                    s = kmer_groups[group_id]['aStart']
                    e = kmer_groups[group_id]['aEnd']
                    kmers[pos][base][kIndex] = {'kseq': kmer_groups[group_id]['kmers'][kIndex], 'aStart': s, 'aEnd': e}

    return kmers


def minimize_kmers(selected_kmers, kmer_groups, clade_info, klen, min_kmers=1, max_kmers=100):
    positions = set(selected_kmers.keys())
    # fix root node conserved kmers '0'
    kmer_pos = {}
    for pos in positions:
        num_in_group = len(selected_kmers[pos]['in_group'])
        num_out_group = len(selected_kmers[pos]['out_group'])

        for base in selected_kmers[pos]['out_group']:
            for kIndex in selected_kmers[pos]['out_group'][base]:
                if not kIndex in kmer_pos:
                    kmer_pos[kIndex] = []
                kmer_pos[kIndex].append(pos)
        for base in selected_kmers[pos]['in_group']:
            for kIndex in selected_kmers[pos]['in_group'][base]:
                if not kIndex in kmer_pos:
                    kmer_pos[kIndex] = []
                kmer_pos[kIndex].append(pos)

    kmer_counts = {}
    for kIndex in kmer_pos:
        kmer_counts[kIndex] = len(kmer_pos[kIndex])
    kmer_counts = {k: v for k, v in sorted(kmer_counts.items(), key=lambda item: item[1])}
    covered_pos = set()
    reduced_kset = set()

    kmer_group_map = {}
    for group_id in kmer_groups:
        for kIndex in kmer_groups[group_id]['kmers']:
            kmer_group_map[kIndex] = group_id

    for kIndex in kmer_counts:
        kmer_id_set = set()
        group_id = kmer_group_map[kIndex]
        for kid in kmer_groups[group_id]['kmers']:
            kmer_id_set.add(kid)
        position_set = set()
        for kid in kmer_id_set:
            position_set = position_set | set(kmer_pos[kid])

        u = position_set - covered_pos
        if len(u) == 0:
            continue
        covered_pos = covered_pos | position_set
        reduced_kset = reduced_kset | kmer_id_set

    for pos in selected_kmers:
        for group in selected_kmers[pos]:
            for base in selected_kmers[pos][group]:
                selected_kmers[pos][group][base] = sorted(list(set(selected_kmers[pos][group][base]) & reduced_kset))




    return selected_kmers


def kmer_worker(outdir, ref_seq, vcf_file, kLen, min_kmer_count, max_kmer_count,prefix, n_threads):
    '''
    Parameters
    ----------
    outdir
    ref_seq
    vcf_file
    kLen
    min_kmer_count
    max_kmer_count
    n_threads
    Returns
    -------
    '''
    logging.info("Creating pseudo-sequences for kmer selection")
    pseudo_seq_file = os.path.join(outdir, "pseudo.seqs.fasta")
    create_pseudoseqs_from_vcf(ref_seq, vcf_file, pseudo_seq_file)

    logging.info("Performing kmer counting")
    # kmer counting
    init_jellyfish_mem = int(os.path.getsize(pseudo_seq_file) / 1024 ** 2 / 10)
    if init_jellyfish_mem == 0:
        init_jellyfish_mem = 1
    jellyfish_file = os.path.join(outdir, "{}-jellyfish.counts.txt".format(prefix))
    run_jellyfish_count(pseudo_seq_file, jellyfish_file, mem="{}M".format(init_jellyfish_mem), k=kLen,
                        n_threads=n_threads)
    if not os.path.isfile(jellyfish_file):
        logging.error("Jellyfish count file failed to generate a result, check logs and try again")
        sys.exit()

    logging.info("Initial kmer filtering")
    kmer_counf_df = parse_jellyfish_counts(jellyfish_file)
    logging.info("Read {} kmers".format(len(kmer_counf_df)))
    kmer_counf_df = kmer_counf_df[kmer_counf_df['count'] >= min_kmer_count]
    kmer_counf_df = kmer_counf_df[kmer_counf_df['count'] <= max_kmer_count]
    kmer_counf_df.reset_index(inplace=True)
    logging.info(" {} kmers remain after count filter >= {}, <= {}".format(len(kmer_counf_df),min_kmer_count,max_kmer_count))
    seqKmers = kmer_counf_df['kmer'].to_dict()
    out_kmer_file = os.path.join(outdir, "{}-filt.kmers.txt".format(prefix))
    logging.info("Collecting k-mer position information")
    return SeqSearchController(seqKmers, pseudo_seq_file, out_kmer_file, n_threads)


def clade_worker(ete_tree_obj, tree_file, vcf_file, metadata_file, min_snp_count, outdir, prefix,max_states=6,
                 min_member_count=1, rcor_thresh=0.4):
    '''
    Parameters
    ----------
    ete_tree_obj
    tree_file
    vcf_file
    min_member_count
    min_snp_count
    max_states
    outdir
    Returns
    -------
    '''
    logging.info("Identifying canonical snps")
    snps = identify_canonical_snps(ete_tree_obj, vcf_file, min_member_count)
    write_snp_report(snps, os.path.join(outdir, "{}-snps.all.txt".format(prefix)))
    snps = snp_based_filter(snps, min_member_count, max_states)
    clade_data = get_valid_nodes(snps, min_snp_count)


    pruned_tree, valid_nodes = prune_tree(ete_tree_obj, list(clade_data.keys()))

    genotypes = generate_genotypes(pruned_tree)
    write_genotypes(genotypes, os.path.join(outdir, "{}-genotypes.supported.txt".format(prefix)), delimeter='.')

    logging.info("Calculating pair-wise distances from nodes to leaves")
    clade_data = calc_node_distances(clade_data, ete_tree_obj)


    logging.info("Parsing sample metadata")
    metadata = parse_metadata(metadata_file)

    if len(metadata) > 0:
        logging.info("Calculating metadata associations")
        clade_data = calc_node_associations(metadata, clade_data, ete_tree_obj)

        if 'year' in metadata[list(metadata.keys())[0]]:
            logging.info("Calculating clade temporal signals")
            clade_data = temporal_signal(metadata, clade_data, ete_tree_obj, rcor_thresh)

        write_node_report(clade_data, os.path.join(outdir, "{}-clades.info.txt".format(prefix)))

    logging.info("Summarizing sample distances")
    matrix_file = os.path.join(outdir, "samples.dists.matrix.csv")
    tree_to_distance_matrix(tree_file, matrix_file)
    sample_dists = get_pairwise_distances_from_matrix(matrix_file)
    max_bins = 100
    dist_summary = sample_dist_summary(sample_dists, num_bins=max_bins)
    nomenclature_dist_ranges = [
        (dist_summary['summary']['dist']['min'], dist_summary['summary']['dist']['25%'] / 2),
        (dist_summary['summary']['dist']['25%'] / 2, dist_summary['summary']['dist']['25%']),
        (dist_summary['summary']['dist']['25%'], dist_summary['summary']['dist']['25%'] +
         dist_summary['summary']['dist']['50%'] / 2),
        (dist_summary['summary']['dist']['25%'] +
         dist_summary['summary']['dist']['50%'] / 2, dist_summary['summary']['dist']['50%']),
        (dist_summary['summary']['dist']['50%'],
         dist_summary['summary']['dist']['50%'] + dist_summary['summary']['dist']['75%'] / 2),
        (dist_summary['summary']['dist']['50%'] + dist_summary['summary']['dist']['75%'] / 2,
         dist_summary['summary']['dist']['75%']),
        (dist_summary['summary']['dist']['75%'],
         dist_summary['summary']['dist']['75%'] + dist_summary['summary']['dist']['max'] / 2),
        (dist_summary['summary']['dist']['75%'] + dist_summary['summary']['dist']['max'] / 2,
         dist_summary['summary']['dist']['75%'] + dist_summary['summary']['dist']['max']),
    ]

    genotypes = generate_genotypes(pruned_tree)
    terminal_nodes = {}
    for sample_id in genotypes:
        node_id = genotypes[sample_id][-1]
        if not node_id in terminal_nodes:
            terminal_nodes[node_id] = 0
        terminal_nodes[node_id] += 1


    candidate_nodes = select_nodes(clade_data, nomenclature_dist_ranges)
    bifurcating_nodes =  set(select_bifucating_nodes(pruned_tree, clade_data))
    valid_nodes = sorted(list(create_compressed_hierarchy(ete_tree_obj, candidate_nodes) | bifurcating_nodes))

    pruned_tree, valid_nodes = prune_tree(ete_tree_obj, valid_nodes)
    genotypes = generate_genotypes(pruned_tree)
    for sample_id in genotypes:
        node_id = genotypes[sample_id][-1]
        if not node_id in terminal_nodes:
            terminal_nodes[node_id] = 0
        terminal_nodes[node_id] += 1
    valid_nodes = set('0')
    for node_id in terminal_nodes:
        if terminal_nodes[node_id] >= min_member_count:
            valid_nodes.add(node_id)

    pruned_tree, valid_nodes = prune_tree(ete_tree_obj, valid_nodes)
    genotypes = generate_genotypes(pruned_tree)
    write_genotypes(genotypes, os.path.join(outdir, "{}-genotypes.selected.txt".format(prefix)), delimeter='.')

    leaf_meta = {}
    for sample_id in genotypes:
        genotype = '.'.join(genotypes[sample_id])
        leaf_meta[sample_id] = [genotype]
        for field in metadata[sample_id]:
            leaf_meta[sample_id].append(metadata[sample_id][field])

    interval_labels = []
    for i in dist_summary['results']['interval'].tolist():
        interval_labels.append("[{}, {})".format(i.left, i.right))
    fig = plot_bar(interval_labels, dist_summary['results']['count'].tolist())

    fig.update_layout(xaxis_title="SNP distance (Hamming)", yaxis_title="Count")
    if fig is not None:
        fig.update_layout(xaxis_title="Intra-clade sample distance", yaxis_title="Count")
        fig = fig.to_html()
        fh = open(os.path.join(outdir, "{}-clade.snp.histo.html".format(prefix)), 'w')
        fh.write(fig)
        fh.close()

    logging.info("Summarizing node snp support")
    node_snp_counts = node_list_attr(clade_data, 'pos')

    if len(set(node_snp_counts)) > 1:
        node_snp_summary = sample_dist_summary(node_snp_counts)
        interval_labels = []
        for i in node_snp_summary['results']['interval'].tolist():
            interval_labels.append("[{}, {})".format(i.left, i.right))
        fig = plot_bar(interval_labels, node_snp_summary['results']['count'].tolist())

        if fig is not None:
            fig.update_layout(xaxis_title="SNP distance (Hamming)", yaxis_title="Count")
            fig = fig.to_html()
            fh = open(os.path.join(outdir, "{}-clade.snp.histo.html".format(prefix)), 'w')
            fh.write(fig)
            fh.close()
    for clade_id in clade_data:
        clade_data[clade_id]['is_selected'] = False
        if clade_id in valid_nodes:
            clade_data[clade_id]['is_selected'] = True

    return clade_data


def metadata_worker(metadata_file, clade_data, ete_tree_obj, outdir, prefix,rcor_thresh=0.2):
    '''
    Parameters
    ----------
    metadata_file
    clade_data
    ete_tree_obj
    outdir
    rcor_thresh
    Returns
    -------
    '''
    logging.info("Parsing sample metadata")
    metadata = parse_metadata(metadata_file)

    if len(metadata) > 0:
        logging.info("Calculating metadata associations")
        clade_data = calc_node_associations(metadata, clade_data, ete_tree_obj)

        if 'year' in metadata[list(metadata.keys())[0]]:
            logging.info("Calculating clade temporal signals")
            clade_data = temporal_signal(metadata, clade_data, ete_tree_obj, rcor_thresh)

        write_node_report(clade_data, os.path.join(outdir, "{}-clades.info.txt".format(prefix)))
    return clade_data

def find_overlaping_gene_feature(start,end,ref_info,ref_name):
    cds_start = start
    cds_end = end
    if not ref_name in ref_info:
        return None
    for feat in ref_info[ref_name]['features']['CDS']:
        positions = feat['positions']
        gene_start = -1
        gene_end = -1
        for s, e in positions:
            if gene_start == -1:
                gene_start = s
            if gene_end < e:
                gene_end = e
            if cds_start  >= s and cds_end  <= e:
                return feat
    return None

def create_scheme_obj(header, selected_kmers, clade_info, sample_genotypes, kmer_rule_obj,ref_features={},trans_table=11):
    '''
    Parameters
    ----------
    header
    selected_kmers
    clade_info
    sample_genotypes
    ref_features
    Returns
    -------
    '''

    perf_annotation = True
    ref_id = list(ref_features.keys())[0]
    ref_seq = ''
    if not 'features' in ref_features[ref_id]:
        perf_annotation = False
        ref_seq = ref_features[ref_id]
    else:
        ref_seq = ref_features[ref_id]['features']['source']

    kmer_key = 0
    unique_genotypes = set()
    for sample_id in sample_genotypes:
        genotype = '.'.join(sample_genotypes[sample_id])
        if not genotype in unique_genotypes:
            unique_genotypes.add(genotype)

    node_members = {}
    for genotype in unique_genotypes:
        genotype = genotype.split('.')
        num_genotypes = len(genotype)
        for i in range(0, num_genotypes):
            node_id = genotype[i]
            if not node_id in node_members:
                node_members[node_id] = ['.'.join(genotype)]
            for k in range(i, num_genotypes):
                g = '.'.join(genotype[0:k + 1])
                node_members[node_id].append(g)
    unique_genotypes = set()
    for node_id in node_members:
        unique_genotypes = unique_genotypes | set(node_members[node_id])

    for node_id in node_members:
        node_members[node_id] = sorted(list(set(node_members[node_id])))

    num_genotypes = len(unique_genotypes)
    scheme = []
    max_entropy = -1
    for pos in selected_kmers:
        ref_base = ref_seq[pos]
        bases = list(selected_kmers[pos].keys())
        alt_bases = []
        for b in bases:
            if b == ref_base:
                continue
            alt_bases.append(b)

        if len(alt_bases) == 0:
            alt_bases = [ref_base]
            #Skip positions where there are multiple kmers to represent a conserved kmer
            if len(selected_kmers[pos][ref_base]) > 1:
                continue
        for i in range(0, len(alt_bases)):
            alt_base = alt_bases[i]
            mutation_key = "snp_{}_{}_{}".format(ref_base, pos + 1, alt_base)
            for k in range(0, len(bases)):
                base = bases[k]
                dna_name = "{}{}{}".format(ref_base, pos + 1, base)
                for kIndex in selected_kmers[pos][base]:
                    obj = {}
                    for field_id in header:
                        obj[field_id] = ''

                    start = selected_kmers[pos][base][kIndex]['aStart']
                    end = selected_kmers[pos][base][kIndex]['aEnd']
                    if perf_annotation:
                        gene_feature = find_overlaping_gene_feature(start, end, ref_features, ref_id)
                    else:
                        gene_feature = None
                        gene_name = 'Intergenic'

                    ref_var_aa = ''
                    alt_var_aa = ''
                    is_cds = False
                    gene_name = ''
                    is_silent = True
                    aa_name = ''
                    aa_var_start = -1
                    if gene_feature is not None:
                        is_cds = True
                        gene_name = gene_feature['gene_name']
                        positions = gene_feature['positions']
                        gene_start = -1
                        gene_end = -1
                        for s, e in positions:
                            if gene_start == -1:
                                gene_start = s
                                gene_end = e
                            else:
                                if gene_start > s:
                                    gene_start = s
                                if gene_end < e:
                                    gene_end = e
                        gene_seq = ref_seq[gene_start:gene_end + 1]
                        rel_var_start = abs(start - gene_start)
                        aa_var_start = int(rel_var_start / 3)
                        codon_var_start = aa_var_start * 3
                        codon_var_end = codon_var_start + 3
                        ref_var_dna = gene_seq[codon_var_start:codon_var_end]
                        ref_var_aa = str(Seq(ref_var_dna).translate(table=trans_table))
                        mod_gene_seq = list(ref_seq[gene_start:gene_end + 1])
                        mod_gene_seq[rel_var_start] = base
                        alt_var_dna = ''.join(mod_gene_seq[codon_var_start:codon_var_end])
                        alt_var_aa = str(Seq(alt_var_dna).translate(table=trans_table))
                        aa_name = "{}{}{}".format(ref_var_aa,aa_var_start+1,alt_var_aa)
                        if alt_var_aa != ref_var_aa:
                            is_silent = False

                    obj['key'] = kmer_key
                    obj['mutation_key'] = mutation_key
                    obj['dna_name'] = dna_name
                    obj['variant_start'] = pos + 1
                    obj['variant_end'] = pos + 1
                    obj['kmer_start'] = selected_kmers[pos][base][kIndex]['aStart'] + 1
                    obj['kmer_end'] = selected_kmers[pos][base][kIndex]['aEnd'] + 1
                    obj['target_variant'] = base
                    obj['target_variant_len'] = 1
                    obj['mutation_type'] = 'snp'
                    obj['ref_state'] = ref_base
                    obj['alt_state'] = alt_base
                    state = 'ref'
                    if base == alt_base and base != ref_base:
                        state = 'alt'
                    obj['state'] = state
                    obj['kseq'] = selected_kmers[pos][base][kIndex]['kseq']
                    obj['klen'] = len(obj['kseq'])
                    obj['homopolymer_len'] = calc_homopolymers(obj['kseq'])
                    obj['is_ambig_ok'] = True
                    obj['is_kmer_found'] = True
                    obj['is_kmer_length_ok'] = True
                    obj['is_kmer_unique'] = True
                    obj['is_valid'] = True
                    obj['kmer_entropy'] = -1
                    obj['positive_genotypes'] = ",".join(kmer_rule_obj[obj['kseq']]['positive_genotypes'])
                    obj['partial_genotypes'] = ",".join(kmer_rule_obj[obj['kseq']]['partial_genotypes'])
                    obj['gene_name'] = gene_name
                    if gene_feature is not None:
                        obj['gene_start'] = gene_start+1
                        obj['gene_end'] = gene_end+1
                        obj['cds_start'] = codon_var_start+1
                        obj['cds_end'] = codon_var_end
                        obj['aa_name'] = aa_name
                        obj['aa_start'] = aa_var_start + 1
                        obj['aa_end'] = aa_var_start + 1


                    obj['is_cds'] = is_cds
                    obj['is_frame_shift'] = False
                    obj['is_silent'] = is_silent

                    is_found = False
                    for clade_id in clade_info:
                        if not clade_info[clade_id]['is_selected']:
                            continue
                        for idx, p in enumerate(clade_info[clade_id]['pos']):
                            if pos != p:
                                continue
                            b = clade_info[clade_id]['bases'][idx]
                            if b != base:
                                continue
                            if clade_id in node_members:
                                counts = [0] * num_genotypes
                                for k in range(0, len(node_members[clade_id])):
                                    counts[k] += 1

                                obj['kmer_entropy'] = calc_shanon_entropy(counts)
                                if obj['kmer_entropy'] > max_entropy:
                                    max_entropy = obj['kmer_entropy']
                                is_found = True
                                break
                        if is_found:
                            break
                    scheme.append(obj)
                    kmer_key += 1
    for i in range(0, len(scheme)):
        e = scheme[i]['kmer_entropy']
        if e == -1:
            scheme[i]['kmer_entropy'] = max_entropy

    return scheme

def construct_ruleset_bck(selected_kmers, genotype_map,outdir,prefix,min_perc):
    kmer_rules = {}
    for pos in selected_kmers:
        for base in selected_kmers[pos]:
            for kIndex in selected_kmers[pos][base]:
                kseq = selected_kmers[pos][base][kIndex]['kseq']
                kmer_rules[kseq] = {}

    genotype_counts = {}
    for sample_id in genotype_map:
        genotype = genotype_map[sample_id].split('.')
        query = []
        for g in genotype:
            query.append(g)
            q = ".".join(query)
            if not q in genotype_counts:
                genotype_counts[q] = 0
            genotype_counts[q]+=1


    kmer_file = os.path.join(outdir, "{}-filt.kmers.txt".format(prefix))
    fh = open(kmer_file,'r')
    header = next(fh).split("\t")
    for line in fh:
        line = line.split("\t")
        if len(line) != 5:
            continue
        kseq= line[3]
        samples = line[4].split(',')
        genotypes = {}
        for sample_id in samples:
            sample_id = sample_id.split("~")[0]
            genotype = genotype_map[sample_id].split('.')
            query = []
            for g in genotype:
                query.append(g)
                q = ".".join(query)
                if not q in genotypes:
                    genotypes[q] = 0
                genotypes[q] += 1

        kmer_rules[kseq] = {'positive_genotypes':[],'partial_genotypes':[]}
        for genotype in genotypes:
            perc = genotypes[genotype] / genotype_counts[genotype]
            if perc >= min_perc:
                kmer_rules[kseq]['positive_genotypes'].append(genotype)
            else:
                kmer_rules[kseq]['partial_genotypes'].append(genotype)

    return kmer_rules

def construct_ruleset(selected_kmers, genotype_map,outdir,prefix,min_perc):
    kmer_rules = {}
    for pos in selected_kmers:
        for base in selected_kmers[pos]:
            for kIndex in selected_kmers[pos][base]:
                kseq = selected_kmers[pos][base][kIndex]['kseq']
                kmer_rules[kseq] = {}

    genotype_counts = {}
    for sample_id in genotype_map:
        genotype = genotype_map[sample_id]
        if not genotype in genotype_counts:
            genotype_counts[genotype] = 0
        genotype_counts[genotype] += 1


    kmer_file = os.path.join(outdir, "{}-filt.kmers.txt".format(prefix))
    fh = open(kmer_file,'r')
    header = next(fh).split("\t")
    for line in fh:
        line = line.split("\t")
        if len(line) != 5:
            continue
        kseq= line[3]
        samples = line[4].split(',')
        genotypes = {}
        for sample_id in samples:
            sample_id = sample_id.split("~")[0]
            genotype = genotype_map[sample_id]
            if not genotype in genotypes:
                genotypes[genotype] = 0

            genotypes[genotype] += 1

        kmer_rules[kseq] = {'positive_genotypes':[],'partial_genotypes':[]}
        for genotype in genotypes:
            perc = genotypes[genotype] / genotype_counts[genotype]
            if perc >= min_perc:
                kmer_rules[kseq]['positive_genotypes'].append(genotype)
            else:
                kmer_rules[kseq]['partial_genotypes'].append(genotype)

    return kmer_rules

def call_consensus_snp_genotypes(ref_seq,vcf_file, genotype_map,outfile,min_perc=1):
    variants = get_variants(vcf_file)
    genotype_counts = {}
    genotype_members = {}
    for sample_id in variants:
        genotype = genotype_map[sample_id]
        if not genotype in genotype_counts:
            genotype_counts[genotype] = {}
            genotype_members[genotype] = 0
        genotype_members[genotype]+=1
        for chrom in variants[sample_id]:
            for pos in variants[sample_id][chrom]:
                base =  variants[sample_id][chrom][pos]
                if not pos in genotype_counts[genotype]:
                    genotype_counts[genotype][pos] = {'A':0,'T':0,'C':0,'G':0,'N':0}
                if not base in genotype_counts[genotype][pos]:
                    base = 'N'
                genotype_counts[genotype][pos][base]+=1
    fh = open(outfile,'w')
    chrom = list(ref_seq.keys())[0]

    for genotype in genotype_counts:
        seq = list(ref_seq[chrom])
        for pos in genotype_counts[genotype]:
            total = genotype_members[genotype]
            b = max(genotype_counts[genotype][pos], key=genotype_counts[genotype][pos].get)
            value = genotype_counts[genotype][pos][b]
            if value / total >= min_perc:
                seq[pos-1] = b
            else:
                seq[pos-1] = 'N'
        fh.write(">{}\n{}\n".format(genotype,"".join(seq)))

    return


def run():
    cmd_args = parse_args()
    tree_file = cmd_args.in_nwk
    variant_file = cmd_args.in_var
    metadata_file = cmd_args.in_meta
    reference_file = cmd_args.reference
    prefix = cmd_args.prefix
    outdir = cmd_args.outdir
    root_method = cmd_args.root_method
    root_name = cmd_args.root_name
    klen = cmd_args.klen
    rcor_thresh = cmd_args.rcor_thresh
    min_snp_count = cmd_args.min_snp_count
    min_member_count = cmd_args.min_members
    min_perc = cmd_args.min_perc
    max_states = cmd_args.max_states
    num_threads = cmd_args.num_threads
    keep_tmp = cmd_args.keep_tmp

    logging = init_console_logger(3)
    os.environ['RAY_worker_register_timeout_seconds'] = '60'
    num_cpus = psutil.cpu_count(logical=False)
    if num_threads > num_cpus:
        num_threads = num_cpus
    tmp_dir_name = tempfile.TemporaryDirectory().name
    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True, num_cpus=num_threads)
    if '.gbk' in reference_file or '.gb' in reference_file:
        seq_file_type = 'genbank'
    else:
        seq_file_type = 'fasta'

    files = [tree_file, variant_file, metadata_file, reference_file]
    for file in files:
        status = validate_file(file)
        if status == False:
            logging.error("Error file {} either does not exist or is empty".format(file))
            sys.exit()

    if root_method is None and root_name is not None:
        root_method = 'outgroup'

    # validate samples present in all files match
    ete_tree_obj = parse_tree(tree_file, logging, ete_format=1, set_root=True, resolve_polytomy=True, ladderize=True,
                              method='midpoint')
    tree_samples = set(ete_tree_obj.get_leaf_names())
    if root_method == 'midpoint':
        root_name = ete_tree_obj.get_tree_root().name

    if root_name not in tree_samples:
        logging.error("Error specified root {} is not in tree: {}".format(root_name, tree_samples))
        sys.exit()

    if root_method == 'outgroup':
        ete_tree_obj = parse_tree(tree_file, logging, ete_format=1, set_root=True, resolve_polytomy=True,
                                  ladderize=True,
                                  method='outgroup',outgroup=root_name)


    for sample_id in tree_samples:
        if isint(sample_id):
            logging.error("Error sample_ids cannot be integers offending sample '{}'".format(sample_id))
            sys.exit()

    vcf_samples = set(vcfReader(variant_file).samples)

    metadata = parse_metadata(metadata_file)

    sample_set = tree_samples | vcf_samples | set(metadata.keys())
    num_samples = len(sample_set)
    missing_samples = sample_set - tree_samples
    if len(missing_samples) > 0:
        logging.error(
            "Error {} samples are not present in tree file: {}".format(len(missing_samples), ','.join(missing_samples)))
        logging.error(
            "Missing samples present in Metadata file {} : {}".format(metadata_file,
                                                                      ','.join(set(metadata.keys()) & missing_samples)))
        logging.error(
            "Missing samples present in Variant file {} : {}".format(variant_file,
                                                                     ','.join(vcf_samples & missing_samples)))
        sys.exit()

    missing_samples = sample_set - vcf_samples
    if len(missing_samples) > 0:
        logging.error("Error {} samples are not present in variant file: {}".format(len(missing_samples),
                                                                                    ','.join(missing_samples)))
        sys.exit()
    missing_samples = sample_set - set(metadata.keys())
    if len(missing_samples) > 0:
        logging.error("Error {} samples are not present in metadata file: {}".format(len(missing_samples),
                                                                                     ','.join(missing_samples)))
        sys.exit()

    if min_member_count > len(sample_set):
        logging.error("Error specified minimum memers to {} but there are only {} samples".format(min_member_count,
                                                                                                  len(sample_set)))
        sys.exit()


    ref_seq = {}
    ref_features = {}
    if seq_file_type == 'genbank':
        ref_features = parse_reference_gbk(reference_file)
        for chrom in ref_features:
            ref_seq[chrom] = ref_features[chrom]['features']['source']
    else:
        with open(reference_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                id = str(record.id)
                seq = str(record.seq).upper()
                ref_seq[id] = seq
        handle.close()

    if not os.path.isdir(outdir):
        logging.info("Creating analysis directory {}".format(outdir))
        os.mkdir(outdir, 0o755)

    analysis_dir = os.path.join(outdir, "_tmp")
    if not os.path.isdir(analysis_dir):
        logging.info("Creating temporary analysis directory {}".format(analysis_dir))
        os.mkdir(analysis_dir, 0o755)

    genotypes = generate_genotypes(ete_tree_obj)
    write_genotypes(genotypes, os.path.join(outdir, "{}-genotypes.raw.txt".format(prefix)), delimeter='.')

    clade_data = clade_worker(ete_tree_obj, tree_file, variant_file, metadata_file, min_snp_count, outdir,prefix , max_states,
                              min_member_count, rcor_thresh)

    valid_nodes = set()
    for clade_id in clade_data:
        if clade_data[clade_id]['is_selected']:
            valid_nodes.add(clade_id)

    pruned_tree, valid_nodes = prune_tree(ete_tree_obj, valid_nodes)
    genotypes = generate_genotypes(pruned_tree)

    leaf_meta = {}
    for sample_id in genotypes:
        genotype = '.'.join(genotypes[sample_id])
        leaf_meta[sample_id] = [genotype]
        for field in metadata[sample_id]:
            leaf_meta[sample_id].append(metadata[sample_id][field])

    logging.info("Creating full tree figure")
    m_ete_tree_obj = copy.deepcopy(ete_tree_obj)
    tree_fig = annotate_tree(os.path.join(outdir, "{}-fulltree.pdf".format(prefix)), ete_tree_obj, valid_nodes, leaf_meta=leaf_meta,h=6400,w=6400,dpi=10000)
    logging.info("Creating representative tree figure")
    plot_single_rep_tree(os.path.join(outdir, "{}-reducedtree.pdf".format(prefix)),  m_ete_tree_obj, valid_nodes, leaf_meta=leaf_meta,h=6400,w=6400,dpi=10000)
    logging.info("Identifying genotyping kmers")
    kmer_groups = kmer_worker(outdir, ref_seq, variant_file, klen, min_member_count, num_samples, prefix, num_threads)

    selected_kmers = select_kmers(kmer_groups, clade_data, klen)


    valid_postions = set(selected_kmers.keys())
    # Preserve the root
    filt = {'0': clade_data['0']}


    for clade_id in clade_data:
        positions = clade_data[clade_id]['pos']
        if len(set(positions) & valid_postions) == 0:
            continue
        chr = []
        fb = []
        pos = []
        for idx, p in enumerate(positions):
            if p in valid_postions:
                chr.append(clade_data[clade_id]['chr'][idx])
                fb.append(clade_data[clade_id]['bases'][idx])
                pos.append(p)
        clade_data[clade_id]['chrom'] = chr
        clade_data[clade_id]['pos'] = pos
        clade_data[clade_id]['bases'] = fb

        filt[clade_id] = clade_data[clade_id]
    clade_data = filt
    logging.info("Readinf genotype assigments")
    logging.info("Readinf genotype assigments")
    genotype_assignments = pd.read_csv(os.path.join(outdir, "{}-genotypes.selected.txt".format(prefix)),sep="\t",header=0).astype(str)
    genotype_assignments = dict(zip(genotype_assignments['sample_id'], genotype_assignments['genotype']))
    logging.info("Constructing kmer rule set")
    kmer_rule_obj = construct_ruleset(selected_kmers, genotype_assignments,outdir,prefix,min_perc)

    logging.info("Creating scheme")
    if len(ref_features) > 0:
        scheme = create_scheme_obj(SCHEME_HEADER, selected_kmers, clade_data, genotypes, kmer_rule_obj,ref_features,)
    else:
        scheme = create_scheme_obj(SCHEME_HEADER, selected_kmers, clade_data, genotypes, kmer_rule_obj,ref_seq)

    call_consensus_snp_genotypes(ref_seq, variant_file, genotype_assignments, os.path.join(outdir,"{}-genotype.consenus.fasta".format(prefix)), min_perc)

    fh = open(os.path.join(outdir, "{}-scheme.txt".format(prefix)), 'w')
    fh.write("{}\n".format("\t".join(SCHEME_HEADER)))
    for i in range(0, len(scheme)):
        row = "\t".join([str(x) for x in list(scheme[i].values())])
        fh.write("{}\n".format(row))
    fh.close()


    if not keep_tmp:
        logging.info("Removing temporary analysis folder")
        shutil.rmtree(analysis_dir)

    logging.info("Analysis complete")