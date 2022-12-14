import os
import shutil
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)

import pandas as pd
import psutil
import ray
from Bio import SeqIO
from Bio.Seq import Seq

from cladeomatic.clades import clade_worker
from cladeomatic.constants import SCHEME_HEADER
from cladeomatic.kmers import kmer_worker
from cladeomatic.utils import init_console_logger, calc_shanon_entropy
from cladeomatic.utils import parse_metadata
from cladeomatic.utils.phylo_tree import parse_tree
from cladeomatic.utils.seqdata import create_pseudoseqs_from_vcf, parse_reference_gbk, calc_homopolymers
from cladeomatic.utils.snpdists import run_snpdists
from cladeomatic.utils.vcfhelper import vcfReader
from cladeomatic.utils.visualization import create_dist_histo
from cladeomatic.version import __version__
from cladeomatic.writers import write_snp_report, write_genotypes, write_node_report, print_params, write_scheme


def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Clade-O-Matic: Genotyping scheme development v. {}".format(__version__),
        formatter_class=CustomFormatter)
    parser.add_argument('--in_var', type=str, required=True,
                        help='Either Variant Call SNP data (.vcf) or TSV SNP data (.txt)')
    parser.add_argument('--in_nwk', type=str, required=False, help='Newick Tree of strains')
    parser.add_argument('--in_groups', type=str, required=False, help='Tab delimited file of genotypes')
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
    parser.add_argument('--max_snp_count', type=int, required=False,
                        help='Maximum number of snps to be selected for each genotype',
                        default=-1)
    parser.add_argument('--min_perc', type=float, required=False,
                        help='Minimum percentage of clade members to be positive for a kmer to be valid', default=1)
    parser.add_argument('--max_states', type=int, required=False,
                        help='Maximum number of states for a position [A,T,C,G,N,-]',
                        default=6)
    parser.add_argument('--max_ambig', type=int, required=False,
                        help='Maximum number of ambiguous bases allowed in a kmer',
                        default=0)
    parser.add_argument('--rcor_thresh', type=float, required=False, help='Correlation coefficient threshold',
                        default=0.4)
    parser.add_argument('--delim', type=str, required=False, help='Genotype delimiter in group file',
                        default=None)
    parser.add_argument('--num_threads', type=int, required=False, help='Number of threads to use', default=1)
    parser.add_argument('--no_plots', required=False, help='Disable plotting', action='store_true')
    parser.add_argument('--keep_tmp', required=False, help='Keep interim files', action='store_true')
    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')
    parser.add_argument('--resume', required=False, help='Resume previous analysis', action='store_true')
    parser.add_argument('--no_compression', required=False, help='Skip compression of tree heirarchy', action='store_true')
    parser.add_argument('--force', required=False, help='Force overwite of existing results directory',
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()

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

def get_tree_genotypes(tree):
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

def parse_group_file(file, delim=None):
    '''

    Parameters
    ----------
    file: TSV file which contains two columns [sample_id, genotype]
    delim: character to split genotypes into levels

    Returns dict of all of the clade memberships in the file
    -------

    '''
    df = pd.read_csv(file, sep="\t", header=0)
    df = df.astype(str)
    groups = {}
    genotypes = df['genotype'].tolist()

    # Determine delimiter
    if delim == None:
        delim_counts = {
            '-': 0,
            ':': 0,
            ';': 0,
            '.': 0,
            '/': 0,
            '|': 0,
            ',': 0,
            ' ': 0,
        }
        for genotype in genotypes:
            for d in delim_counts:
                delim_counts[d] += genotype.count(d)
        delim = max(delim_counts, key=delim_counts.get)

    valid_nodes = set()
    # Build genotype lookup
    for genotype in genotypes:
        genotype = genotype.split(delim)
        for node_id in [str(x) for x in genotype]:
            valid_nodes.add(node_id)
        for i in range(1, len(genotype) + 1):
            g = "{}".format(delim).join([str(x) for x in genotype[0:i]])
            groups[g] = set()

    sample_map = {}
    # Add sample to each genotype
    for row in df.itertuples():
        sample_id = row.sample_id
        sample_map[row.index] = {'sample_id': sample_id, 'genotype': row.genotype}
        genotype = row.genotype.split(delim)
        for i in range(1, len(genotype) + 1):
            g = "{}".format(delim).join([str(x) for x in genotype[0:i]])
            groups[g].add(row.index)

    return {'delimeter': delim, 'sample_map': sample_map, 'membership': groups,
            'genotypes': set(df['genotype'].tolist()), 'valid_nodes': valid_nodes}

def parse_tree_groups(ete_tree_obj, delim='.'):
    '''

    Parameters
    ----------
    ete_tree_obj: ETE3 tree object
    delim: character to split genotypes into levels

    Returns dict of all of the clade memberships in the file
    -------

    '''
    genotypes = get_tree_genotypes(ete_tree_obj)

    id = 0
    sample_map = {}
    groups = {}
    unique_genotypes = set()
    valid_nodes = set()
    for sample_id in genotypes:
        genotype = "{}".format(delim).join(genotypes[sample_id])
        unique_genotypes.add(genotype)
        sample_map[id] = {'sample_id': sample_id, 'genotype': genotype}
        genotype = genotype.split(delim)
        for node_id in [str(x) for x in genotype]:
            valid_nodes.add(node_id)
        for i in range(1, len(genotype) + 1):
            g = "{}".format(delim).join([str(x) for x in genotype[0:i]])
            groups[g] = set()
        id += 1

    # Add sample to each genotype
    for id in sample_map:
        genotype = sample_map[id]['genotype'].split(delim)
        for i in range(1, len(genotype) + 1):
            g = "{}".format(delim).join([str(x) for x in genotype[0:i]])
            groups[g].add(id)

    return {'delimeter': delim, 'sample_map': sample_map, 'membership': groups,
            'genotypes': unique_genotypes, 'valid_nodes': valid_nodes}

def find_overlaping_gene_feature(start, end, ref_info, ref_name):
    '''

    Parameters
    ----------
    start
    end
    ref_info
    ref_name

    Returns
    -------

    '''
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
            if cds_start >= s and cds_end <= e:
                return feat
    return None

def create_scheme(header,ref_features,kmer_worker,sample_genotypes,trans_table=11):
    perf_annotation = True
    ref_id = list(ref_features.keys())[0]
    ref_seq = ''
    if not 'features' in ref_features[ref_id]:
        perf_annotation = False
        ref_seq = ref_features[ref_id]
    else:
        ref_seq = ref_features[ref_id]['features']['source']


    selected_kmers = kmer_worker.kmer_scheme_data
    kmer_rule_obj = kmer_worker.rule_set
    kmer_info = kmer_worker.int_base_kmer_lookup

    unique_genotypes = set(sample_genotypes.values())
    num_genotypes = len(unique_genotypes)
    scheme = []
    max_entropy = -1
    mutation_keys = {}
    kmer_key = 0

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
            # Skip positions where there are multiple kmers to represent a conserved kmer
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
                    kseq = kmer_worker.get_kseq_by_index(kIndex)
                    for field_id in header:
                        obj[field_id] = ''
                    start = kmer_info[kIndex]['aln_start']
                    end = kmer_info[kIndex]['aln_end']
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
                        aa_name = "{}{}{}".format(ref_var_aa, aa_var_start + 1, alt_var_aa)
                        if alt_var_aa != ref_var_aa:
                            is_silent = False

                    obj['key'] = kmer_key
                    obj['mutation_key'] = mutation_key
                    obj['dna_name'] = dna_name
                    obj['variant_start'] = pos + 1
                    obj['variant_end'] = pos + 1
                    obj['kmer_start'] = start + 1
                    obj['kmer_end'] = end + 1
                    obj['target_variant'] = base
                    obj['target_variant_len'] = 1
                    obj['mutation_type'] = 'snp'
                    obj['ref_state'] = ref_base
                    obj['alt_state'] = alt_base
                    state = 'ref'
                    if base == alt_base and base != ref_base:
                        state = 'alt'
                    obj['state'] = state
                    obj['kseq'] = kseq
                    obj['klen'] = len(obj['kseq'])
                    obj['homopolymer_len'] = calc_homopolymers(obj['kseq'])
                    obj['is_ambig_ok'] = True
                    obj['is_kmer_found'] = True
                    obj['is_kmer_length_ok'] = True
                    obj['is_kmer_unique'] = True
                    obj['is_valid'] = True

                    counts = [0] * num_genotypes
                    geno_found = unique_genotypes & (
                        set(kmer_rule_obj[kIndex]['positive_genotypes']) | set(
                            kmer_rule_obj[kIndex]['partial_genotypes']))
                    num_found = len(geno_found)
                    if num_genotypes < num_found:
                        print("{}\t{}\t{}\t{}".format(pos,kIndex,num_genotypes,num_found ))
                        print(unique_genotypes)
                        print(set(kmer_rule_obj[kIndex]['positive_genotypes']))
                        print(set(kmer_rule_obj[kIndex]['partial_genotypes']))
                    for j in range(0, num_found):
                        counts[j] = 1
                    obj['kmer_entropy'] = calc_shanon_entropy(counts)
                    if max_entropy < obj['kmer_entropy']:
                        max_entropy = obj['kmer_entropy']

                    obj['positive_genotypes'] = ",".join(sorted(kmer_rule_obj[kIndex]['positive_genotypes']))
                    obj['partial_genotypes'] = ",".join(sorted(kmer_rule_obj[kIndex]['partial_genotypes']))
                    obj['gene_name'] = gene_name
                    if not mutation_key in mutation_keys:
                        mutation_keys[mutation_key] = {'alt': [], 'ref': []}

                    if len(obj['positive_genotypes']) > 0:
                        mutation_keys[mutation_key][state].append(obj['positive_genotypes'])

                    if gene_feature is not None:
                        obj['gene_start'] = gene_start + 1
                        obj['gene_end'] = gene_end + 1
                        obj['cds_start'] = codon_var_start + 1
                        obj['cds_end'] = codon_var_end
                        obj['aa_name'] = aa_name
                        obj['aa_start'] = aa_var_start + 1
                        obj['aa_end'] = aa_var_start + 1

                    obj['is_cds'] = is_cds
                    obj['is_frame_shift'] = False
                    obj['is_silent'] = is_silent
                    scheme.append(obj)
                    kmer_key += 1

    mutations_keys_to_remove = []
    for mkey in mutation_keys:
        if len(mutation_keys[mkey]['ref']) == 0 and len(mutation_keys[mkey]['alt']) == 0:
            mutations_keys_to_remove.append(mkey)
    filt = []
    for i in range(0, len(scheme)):
        if len(scheme[i]['positive_genotypes']) == 0 and len(scheme[i]['partial_genotypes']) == 0:
            scheme[i]['kmer_entropy'] = max_entropy
        if scheme[i]['mutation_key'] in mutations_keys_to_remove:
            continue
        filt.append(scheme[i])
    return scheme

def run():
    cmd_args = parse_args()
    tree_file = cmd_args.in_nwk
    group_file = cmd_args.in_groups
    variant_file = cmd_args.in_var
    metadata_file = cmd_args.in_meta
    reference_file = cmd_args.reference
    prefix = cmd_args.prefix
    outdir = cmd_args.outdir
    root_method = cmd_args.root_method
    root_name = cmd_args.root_name
    klen = cmd_args.klen
    rcor_thresh = cmd_args.rcor_thresh
    max_snp_count = cmd_args.max_snp_count
    min_snp_count = cmd_args.min_snp_count
    min_member_count = cmd_args.min_members
    min_perc = cmd_args.min_perc
    max_states = cmd_args.max_states
    num_threads = cmd_args.num_threads
    keep_tmp = cmd_args.keep_tmp
    delim = cmd_args.delim
    max_ambig = cmd_args.max_ambig
    no_compression = cmd_args.no_compression
    force = cmd_args.force

    logging = init_console_logger(3)

    # Initialize Ray components
    os.environ['RAY_worker_register_timeout_seconds'] = '60'
    num_cpus = psutil.cpu_count(logical=False)
    if num_threads > num_cpus:
        num_threads = num_cpus

    if not ray.is_initialized():
        ray.init(num_cpus=num_threads)

    # Initialize directory
    if not os.path.isdir(outdir):
        logging.info("Creating analysis directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    elif not force:
        logging.info("Error directory {} already exists, if you want to overwrite existing results then specify --force".format(outdir))
        sys.exit()

    print_params(cmd_args, os.path.join(outdir,"{}-params.log".format(prefix)))

    analysis_dir = os.path.join(outdir, "_tmp")
    if not os.path.isdir(analysis_dir):
        logging.info("Creating temporary analysis directory {}".format(analysis_dir))
        os.mkdir(analysis_dir, 0o755)

    # Toggle tree based or group based identification
    if tree_file is not None:
        mode = 'tree'
        status = validate_file(tree_file)
        if status == False:
            logging.error("Error file {} either does not exist or is empty".format(tree_file))
            sys.exit()
    elif group_file is not None:
        mode = 'group'
        status = validate_file(group_file)
        if status == False:
            logging.error("Error file {} either does not exist or is empty".format(group_file))
            sys.exit()
    else:
        logging.error("You must specify either a tree file or group file")
        sys.exit()

    # Validate input file presence
    files = [variant_file, metadata_file, reference_file]
    for file in files:
        status = validate_file(file)
        if status == False:
            logging.error("Error file {} either does not exist or is empty".format(file))
            sys.exit()

    # Parse reference sequence
    if '.gbk' in reference_file or '.gb' in reference_file:
        seq_file_type = 'genbank'
    else:
        seq_file_type = 'fasta'

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
    ref_seq_id = list(ref_seq.keys())[0]
    vcf_samples = set(vcfReader(variant_file).samples)
    metadata = parse_metadata(metadata_file)

    if mode == 'tree':
        if root_method is None and root_name is not None:
            root_method = 'outgroup'

        # validate samples present in all files match
        ete_tree_obj = parse_tree(tree_file, logging, ete_format=1, set_root=True, resolve_polytomy=True,
                                  ladderize=True,
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
                                      method='outgroup', outgroup=root_name)

        for sample_id in tree_samples:
            if isint(sample_id):
                logging.error("Error sample_ids cannot be integers offending sample '{}'".format(sample_id))
                sys.exit()
        group_samples = tree_samples
        sample_set = tree_samples | vcf_samples | set(metadata.keys())
        num_samples = len(sample_set)
        missing_samples = sample_set - tree_samples
    else:
        group_data = parse_group_file(group_file, delim=delim)
        group_samples = set(group_data['sample_list'])
        sample_set = group_samples | vcf_samples | set(metadata.keys())
        num_samples = len(sample_set)
        missing_samples = sample_set - group_samples

    # Validate Sample id's accross all inputs
    if len(missing_samples) > 0:
        logging.error(
            "Error {} samples are not present in tree or genotype file: {}".format(len(missing_samples),
                                                                                   ','.join(missing_samples)))
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

    if mode == 'tree':
        group_data = parse_tree_groups(ete_tree_obj)
    else:
        group_data = parse_group_file(group_file)

    logging.info("Reading metadata file")
    metadata = parse_metadata(metadata_file)

    logging.info("Recreating fasta sequences from vcf")
    pseudo_seq_file = os.path.join(outdir, "pseudo.seqs.fasta")
    create_pseudoseqs_from_vcf(ref_seq_id,ref_seq[ref_seq_id], variant_file, pseudo_seq_file)

    #calculate distance matrix
    logging.info("Calculating SNP distance matrix")
    distance_matrix_file = os.path.join(outdir,"{}-dist.mat.txt".format(prefix))
    run_snpdists(pseudo_seq_file, distance_matrix_file, num_threads)

    #perform clade-snp work
    perform_compression = True
    if no_compression:
        perform_compression = False
    logging.info("Performing canonical SNP detection")
    cw = clade_worker(variant_file, metadata , distance_matrix_file, group_data, ref_seq[ref_seq_id], perform_compression=perform_compression,delim=delim,
                      min_snp_count=min_snp_count, max_snps=max_snp_count, max_states=max_states, min_members=min_member_count,
                 min_inter_clade_dist=1, num_threads=num_threads,rcor_thresh=rcor_thresh)

    logging.info("Read {} variant positions from {}".format(cw.num_positions,variant_file))
    logging.info("Found {} valid variant positions".format(cw.num_valid_positions))
    logging.info("Found {} canonical variant positions".format(cw.num_canonical_positions))
    logging.info("Initial set of {} genotyping positions selected".format(len(cw.selected_positions)))

    write_node_report(cw.clade_data, os.path.join(outdir, "{}-clades.info.txt".format(prefix)))
    write_snp_report(cw.snp_data, os.path.join(outdir, "{}-snps.info.txt".format(prefix)))

    write_genotypes(cw.raw_genotypes, os.path.join(outdir, "{}-genotypes.raw.txt".format(prefix)))
    write_genotypes(cw.supported_genotypes, os.path.join(outdir, "{}-genotypes.supported.txt".format(prefix)))


    create_dist_histo(cw.distance_histo,os.path.join(outdir, "{}-sample.distances.html".format(prefix)))


    #perform kmer selection
    genotype_map = cw.selected_genotypes
    target_positions = cw.selected_positions

    logging.info("Performing kmer selection")
    kw = kmer_worker(ref_seq[ref_seq_id], pseudo_seq_file, analysis_dir, prefix, klen, genotype_map, max_ambig=max_ambig, min_perc=min_perc,
                 target_positions=target_positions, num_threads=num_threads)


    #Update based on positions which could not be assigned a kmer
    positions_missing_kmer = kw.positions_missing_kmer
    if len(positions_missing_kmer) > 0:
        cw.remove_snps(list(positions_missing_kmer.keys()))
    logging.info("A total of {} genotyping positions removed due to no valid kmers found".format(len(positions_missing_kmer)))
    logging.info("Final set of {} genotyping positions selected".format(len(cw.selected_positions)))
    write_genotypes(cw.selected_genotypes, os.path.join(outdir, "{}-genotypes.selected.txt".format(prefix)))

    if max_snp_count > 1:
        cw.prune_snps()
        variant_pos = set(cw.variant_positions)
        selected_positions = set(cw.get_selected_positions())
        kw.remove_scheme_pos(variant_pos - selected_positions)

    logging.info("Creating scheme")
    if len(ref_features) > 0:
        scheme = create_scheme(SCHEME_HEADER,ref_features,kw,cw.selected_genotypes,trans_table=11)
    else:
        scheme = create_scheme(SCHEME_HEADER, ref_seq, kw, cw.selected_genotypes, trans_table=11)

    write_scheme(SCHEME_HEADER, scheme, os.path.join(outdir, "{}-scheme.txt".format(prefix)))

    if not keep_tmp:
        logging.info("Removing temporary analysis folder")
        shutil.rmtree(analysis_dir)

    logging.info("Analysis complete")
