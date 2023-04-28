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
from cladeomatic.constants import MIN_FILE_SIZE

def parse_args():
    """ Arugment Parsing method.

    A function to parse the command line arguments passed at initilization of Clade-o-matic,
    format these arguments,  and return help prompts to the user shell when specified.

    Returns
    -------
    ArgumentParser object
        The arguments and their user specifications, the usage help prompts and the correct formatting
        for the incoming argument (str, int, etc.)
    """

    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        """
        Class to instantiate the formatter classes required for the argument parser.
        Required for the correct formatting of the default parser values

        Parameters
        ----------
        ArgumentDefaultsHelpFormatter object
            Instatiates the default values for the ArgumentParser for display on the command line.
        RawDescriptionHelpFormatter object
            Ensures the correct display of the default values for the ArgumentParser
        """

        pass
    #creation of the ArgumentParser object
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
                        help='Minimum number of unique SNPs for a clade to be valid',
                        default=1)
    parser.add_argument('--max_snp_count', type=int, required=False,
                        help='Maximum number of SNPs to be selected for defining each genotype to prevent large numbers of redundant SNPs',
                        default=-1)
    parser.add_argument('--min_perc', type=float, required=False,
                        help='Minimum percentage of clade members to be positive for a kmer to be valid', default=0.1)
    parser.add_argument('--max_site_ambig', type=float, required=False,
                        help='Maximum percentage of input sequences which can be missing a site for it to still be valid', default=0.25)
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
    parser.add_argument('--no_compression', required=False, help='Skip compression of tree hierarchy', action='store_true')
    parser.add_argument('--force', required=False, help='Force overwrite of existing results directory',
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()

def validate_file(file):
    """ File validation method.
    A method to determine if the file path leads to a valid file.  Return true if
    the file exists and is not empty.

    Parameters
    ----------
    file : str
        The path to the file

    Returns
    -------
    bool
        True if the file exists and is not empty, False otherwise
    """

    if os.path.isfile(file) and os.path.getsize(file) > 32:
        return True
    else:
        return False

def isint(str):
    """
    A function to determine if a string can be cast to an integer

    Parameters
    ----------
    str: str
        The string to be tested for integer casting

    Returns
    -------
    bool
        True if the string can be cast to an integer, False if not
    """

    try:
        int(str)
        return True
    except ValueError:
        return False

def get_tree_genotypes(tree):
    """
    This method accepts an ETE3 tree object and returns the
    heirarchical node ids for every leaf.

    Parameters
    ----------
    tree : ETE3 tree object
        The tree to be processed for its genotypes

    Returns
    -------
    dict
        A dictionary of the leaf names and node heirarchy - [node id, list(member ids)]

    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation
    """

    samples = tree.get_leaf_names()
    geno = {}
    for sample_id in samples:
        geno[sample_id] = []
    #Preorder traversal sequence is root, left, right - see more at ete3 documentation
    for node in tree.traverse("preorder"):
        node_id = node.name
        if node.is_leaf():
            continue
        leaves = node.get_leaf_names()
        for sample_id in leaves:
            geno[sample_id].append(node_id)

    return geno

def parse_group_file(file, delim=None):
    """
    Method to read the grouping file passed and construct an ETE3 object from
    values of the grouping file.

    Parameters
    ----------
    file : str
        The path to TSV file which contains two columns [sample_id, genotype]
    delim : char
        The character to split genotypes into levels, default is None as it will be determined dynamically

    Returns
    -------
    dict
        A dictionary of all the clade memberships in the file: The delimeter used, the sample map, the group membership, the genotypes, and valid nodes.

    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation
    """

    df = pd.read_csv(file, sep="\t", header=0)
    df = df.astype(str)
    groups = {}
    genotypes = df['genotype'].tolist()

    # Determine delimiter from the passed data file
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
    # Add sample to each genotype determined
    for row in df.itertuples(index=True):
        sample_id = row.sample_id
        sample_map[sample_id] = {'sample_id': sample_id, 'genotype': row.genotype}
        genotype = row.genotype.split(delim)
        for i in range(1, len(genotype) + 1):
            g = "{}".format(delim).join([str(x) for x in genotype[0:i]])
            groups[g].add(sample_id)

    return {'delimeter': delim, 'sample_map': sample_map, 'membership': groups,
            'genotypes': set(df['genotype'].tolist()), 'valid_nodes': valid_nodes}

def parse_tree_groups(ete_tree_obj, delim='.'):
    """
    Method takes an ETE3 tree object, parses it, constructs the group data, and returns it in the
    form of a dictionary.

    Parameters
    ----------
    ete_tree_obj : ETE3 tree object
        The ETE3 tree object for parsing
    delim : char
        The character to split genotypes into levels, default is '.', note delimiter is contained in the ete3 object.

    Returns
    -------
    dict
        A dictionary of all the clade memberships in the file. The dictionary consists the delimeter used, the sample map,the group membership, the unique genotypes, and the valid nodes.

    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation
    """
    genotypes = get_tree_genotypes(ete_tree_obj)

    id = 0
    sample_map = {}
    groups = {}
    unique_genotypes = set()
    valid_nodes = set()
    #parse out the genotype ids and sample ids from the tree
    for sample_id in genotypes:
        genotype = "{}".format(delim).join(genotypes[sample_id])
        unique_genotypes.add(genotype)
        sample_map[id] = {'sample_id': sample_id, 'genotype': genotype}
        genotype = genotype.split(delim)
        #parse the tree to find the groups
        for node_id in [str(x) for x in genotype]:
            valid_nodes.add(node_id)
        for i in range(1, len(genotype) + 1):
            g = "{}".format(delim).join([str(x) for x in genotype[0:i]])
            groups[g] = set()
        id += 1

    # Add sample to each genotype node
    for id in sample_map:
        genotype = sample_map[id]['genotype'].split(delim)
        for i in range(1, len(genotype) + 1):
            g = "{}".format(delim).join([str(x) for x in genotype[0:i]])
            groups[g].add(id)

    return {'delimeter': delim, 'sample_map': sample_map, 'membership': groups,
            'genotypes': unique_genotypes, 'valid_nodes': valid_nodes}

def find_overlaping_gene_feature(start, end, ref_info, ref_name):
    """
    Method to find the overlapping gene features from the reference
    sequence dictionary passed to the method.

    Parameters
    ----------
    start : int
        The start of the coding sequence
    end : int
        The end of the coding sequence
    ref_info : dict
        The dictionary object containing the reference sequence and associated data (sequences, features, coding regions, and their positions).
    ref_name : str
        The name for the reference sequence.  This dictionary has the gene name, sequence for the gene, the corresponding amino acid sequence, the position, and the gene length.

    Returns
    -------
    dict
        The dictionary for the gene features

    """
    cds_start = start
    cds_end = end
    #if the refrence sequence name is not in the sequence dictionary
    #return None - invalid
    if not ref_name in ref_info:
        return None
    #loop through the reference dictionary
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
    """  Creates the kmer-based scheme.

    This method creates the kmer-based scheme by iterating through
    the kmer worker object to obtain the kmer set, scheme data, rule set,
    genotype list, and reference features (if present).  This method then uses these
    items to construct a kmer scheme containing the SNP and kmer identifiers
    for downstream use.  Please refer to the file
    examples/small_test/cladeomatic/cladeomatic-kmer.scheme.txt for more information.

    Parameters
    ----------
    header : str
        The kmer scheme header row.  The default string is the :const:`cladeomatic.constants.SNP_HEADER`.
    ref_features : dict
        The dictionary containing the reference features
    kmer_worker : kmer_worker object
        The kmer_worker object (:class:`cladeomatic.kmers.kmer_worker`) that contains the kmer set, scheme data, and rule set for kmer scheme creation (kmer_scheme_data, rule_set, int_base_kmer_lookup, get_kseq_by_index)
    sample_genotypes : dict
        The dictionary object that contains the selected genotype identifiers for the samples for use in the scheme
    trans_table : int
        The signifier for which NCBI amino acid translation table you wish to use - default is table 11 for Bacterial, Archaeal and Plant Plastid Code.  Refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more details

    Returns
    -------
    list
        The list for use for the kmer scheme containing all elements of the kmer scheme.  Refer to the sample scheme above.

    """
    perf_annotation = True
    ref_id = list(ref_features.keys())[0]
    ref_seq = ''
    #add feature annotation if specified
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
    #initialization for the maximum value for the shannon's entropy calculation
    max_entropy = -1
    mutation_keys = {}
    kmer_key = 0

    #loop through the selected kmers dictionary
    for pos in selected_kmers:
        ref_base = ref_seq[pos]
        #ensure the reference base for the position in the ref sequence
        #is not in the same position for the kmer - if so skip to avoid
        #kmer duplication
        if not ref_base in selected_kmers[pos]:
            continue
        #get the list of base keys for the kmers identified for the samples
        #submitted
        bases = list(selected_kmers[pos].keys())
        alt_bases = []

        #loop through the list of kmer keys - if the base is the same
        #as the reference base, skip to avoid duplication
        for b in bases:
            if b == ref_base:
                continue
            alt_bases.append(b)

        if len(alt_bases) == 0:
            alt_bases = [ref_base]
            # Skip positions where there are multiple kmers to represent a conserved kmer
            if len(selected_kmers[pos][ref_base]) > 1:
                continue
        #loop through the base keys to retreive details for the scheme
        for i in range(0, len(alt_bases)):
            alt_base = alt_bases[i]
            #format the SNP name for the scheme
            mutation_key = "snp_{}_{}_{}".format(ref_base, pos + 1, alt_base)
            #loop through the base keys to find and append the kmers
            for k in range(0, len(bases)):
                base = bases[k]
                #format the DNA identifier for the scheme
                dna_name = "{}{}{}".format(ref_base, pos + 1, base)
                for kIndex in selected_kmers[pos][base]:
                    #object dictionary for the scheme aspects
                    obj = {}
                    #get the kmer sequence for the identifier passed
                    kseq = kmer_worker.get_kseq_by_index(kIndex)
                    #create the scheme header in the object dictionary
                    for field_id in header:
                        obj[field_id] = ''
                    start = kmer_info[kIndex]['aln_start']
                    end = kmer_info[kIndex]['aln_end']
                    #if there is annotation to be included, retreive it here
                    if perf_annotation:
                        gene_feature = find_overlaping_gene_feature(start, end, ref_features, ref_id)
                    else:
                        gene_feature = None
                        gene_name = 'Intergenic'
                    #initialize the variables for the fields required for the scheme
                    ref_var_aa = ''
                    alt_var_aa = ''
                    is_cds = False
                    gene_name = ''
                    is_silent = True
                    aa_name = ''
                    aa_var_start = -1
                    #add and format the annotations for scheme additions
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
                    #add the variable outputs for the scheme to the object dictionary
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

                    kmer_rule_obj[kIndex]['positive_genotypes'] = list(unique_genotypes & set(kmer_rule_obj[kIndex]['positive_genotypes']))
                    kmer_rule_obj[kIndex]['partial_genotypes'] = list(
                        unique_genotypes & set(kmer_rule_obj[kIndex]['partial_genotypes']))
                    counts = [0] * num_genotypes
                    geno_found =  set(kmer_rule_obj[kIndex]['positive_genotypes']) | set(
                            kmer_rule_obj[kIndex]['partial_genotypes'])
                    num_found = len(geno_found)

                    for j in range(0, num_found):
                        counts[j] = 1
                    #calculate the shannon entropy
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
    #clean up: remove any mutations that may have gotten through
    #the original filtering
    mutations_keys_to_remove = []
    for mkey in mutation_keys:
        if len(mutation_keys[mkey]['ref']) == 0 and len(mutation_keys[mkey]['alt']) == 0:
            mutations_keys_to_remove.append(mkey)

    filt = []
    kmer_key = 0
    for i in range(0, len(scheme)):
        if len(scheme[i]['positive_genotypes']) == 0 and len(scheme[i]['partial_genotypes']) == 0:
            scheme[i]['kmer_entropy'] = max_entropy
        if scheme[i]['mutation_key'] in mutations_keys_to_remove:
            continue
        scheme[i]['key'] = kmer_key
        filt.append(scheme[i])
        kmer_key+=1
    scheme = filt
    return scheme

def filter_vcf(input_vcf,output_vcf,max_states,max_missing):
    """
    This method filters the VCF input file for invalid SNP sites, writes
    a filtered VCF file to the path specified, and returns the number of
    SNPs removed in the process.

    Parameters
    ----------
    input_vcf : str
        The file path to the VCF file to be read for values
    output_vcf : str
        The file path to the filtered VCF file to be written
    max_states : int
        The maximum number of variant states allowed for this analysis
    max_missing : int
        The maximum number of missing variant states allowed for this analysis

    Returns
    -------
    int
        The number of mising SNPs removed from the variant call file.  These SNPs are removed if they violate the processing rules for the VCF file.
    """
    #read the vcf file and return the vcf object: NOTE files
    vcf = vcfReader(input_vcf)
    #get the data from the VCF file in the form of a dictionary
    data = vcf.process_row()
    #retrieve the sample information from the VCF file in the form of a
    #list
    samples = vcf.samples
    num_samples = len(samples)
    if data is None:
        shutil.copy(input_vcf,output_vcf)
        return

    invalid_positions = {}

    while data is not None:
        chrom = data['#CHROM']
        #initialize the invalid positions dictionary for the
        #removal of invalid samples based on bad base calls
        if not chrom in invalid_positions:
            invalid_positions[chrom] = []

        pos = int(data['POS']) - 1
        #initialize the counts for the bases
        base_counts = {
            'A':0,
            'T':0,
            'C':0,
            'G':0,
            'N':0
        }
        #find the bad base calls
        for sample_id in samples:
            base = data[sample_id]
            if base not in ['A','T','C','G']:
                base = 'N'
            base_counts[base]+=1
        #if the ratio of bad calls to the numver of samples is
        #greater than the max allowed missing samples
        #add them to the invalid positions dictionary
        if base_counts['N'] / num_samples > max_missing:
            invalid_positions[chrom].append(pos)

        #count the number of states in the VCF
        count_states = 0
        for base in ['A','T','C','G']:
            if base_counts[base] > 0:
                count_states+=1

        if count_states == 1 or count_states > max_states:
            invalid_positions[chrom].append(pos)
        #because this is python, re-initialize just in case
        data = vcf.process_row()
    #nullify the vcf object instance
    del(vcf)
    #open the input and output files, format and write the filtered output
    in_fh = open(input_vcf,'r')
    out_fh = open(output_vcf,'w')
    count_snps_removed = 0
    for line in in_fh:
        line = line.rstrip()
        row = line.split("\t")
        if len(row) == 0:
            continue
        if len(row) < 2:
            out_fh.write("{}\n".format(line))
            continue
        chrom = row[0]
        if chrom == 'CHROM' or not row[1].isnumeric():
            out_fh.write("{}\n".format(line))
            continue

        pos = int(row[1]) - 1
        if chrom in invalid_positions:
            if pos in invalid_positions[chrom]:
                count_snps_removed+=1
                continue
        out_fh.write("{}\n".format(line))
    in_fh.close()
    out_fh.close()

    return count_snps_removed

def create_snp_scheme(header,ref_features,clade_obj,trans_table=11):
    """ Creates the SNP scheme for file writing

    This method is similar to :meth:`create_scheme`, but it does not include
    the kmer sequences in the SNP scheme. The method iterates through
    the genotype list from the clade object and reference features
    (if present) to construct a scheme containing the SNP identifiers.  Please refer to
    the file examples/small_test/cladeomatic/cladeomatic-snp.scheme.txt for more information.

    Parameters
    ---------
    header : str
        The snp scheme header row.  The default string is the :const:`cladeomatic.constants.SNP_HEADER`.
    ref_features : dict
        The dictionary containing the reference features
    clade_obj : clade_worker object
        The clade_worker object (:class:`cladeomatic.clades.clade_worker`) that contains the clade data, selected genotypes, scheme data, and rule set for snp scheme creation
    sample_genotypes : dict
        The dictionary object that contains the selected genotype identifiers for the samples for use in the scheme
    trans_table : int
        The signifier for which NCBI amino acid translation table you wish to use - default is table 11 for Bacterial, Archaeal and Plant Plastid Code.  Refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more details

    Returns
    -------
    list
        The list for use for the snp scheme containing all elements of the snp scheme.  Refer to the sample scheme above.

    """
    #retrieve the sample genotypes from the cladeworker object
    sample_genotypes = clade_obj.selected_genotypes
    perf_annotation = True
    #retrieve the list of feature keys if they exist
    ref_id = list(ref_features.keys())[0]
    ref_seq = ''
    if not 'features' in ref_features[ref_id]:
        perf_annotation = False
        ref_seq = ref_features[ref_id]
    else:
        ref_seq = ref_features[ref_id]['features']['source']
    #retrieve the set of genotypes from the sample values
    unique_genotypes = set(sample_genotypes.values())
    num_genotypes = len(unique_genotypes)
    #initialize variables
    scheme = []
    max_entropy = -1
    mutation_keys = {}
    row_index = 0
    #retreive the rules
    rules = clade_obj.get_genotype_snp_rules()
    selected_pos_ref_states = clade_obj.selected_pos_ref_states
    valid_bases = ["A","T","C","G"]
    scheme = []
    #loop through the reference dictionary for the scheme data
    for pos in selected_pos_ref_states:
        ref_base = selected_pos_ref_states[pos]
        alt_bases = []
        bases_present = []
        #find the alternate bases aka mutations for scheme output
        for base in valid_bases:
            if len(rules[pos][base]['positive_genotypes'] ) > 0 or len(rules[pos][base]['partial_genotypes'] ):
                bases_present.append(base)
                if base != ref_base:
                    alt_bases.append(base)
        #loop thorugh the alternate bases/mutations to create the scheme output
        for i in range(0, len(alt_bases)):
            alt_base = alt_bases[i]
            mutation_key = "snp_{}_{}_{}".format(ref_base, pos + 1, alt_base)
            for k in range(0, len(bases_present)):
                base = bases_present[k]
                dna_name = "{}{}{}".format(ref_base, pos + 1, base)
                row_obj = {}
                for field_id in header:
                    row_obj[field_id] = ''
                if perf_annotation:
                    gene_feature = find_overlaping_gene_feature(pos, pos, ref_features, ref_id)
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
                    rel_var_start = abs(pos - gene_start)
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

                row_obj['key'] = row_index
                row_obj['mutation_key'] = mutation_key
                row_obj['dna_name'] = dna_name
                row_obj['variant_start'] = pos + 1
                row_obj['variant_end'] = pos + 1
                row_obj['kmer_start'] = -1
                row_obj['kmer_end'] = -1
                row_obj['target_variant'] = base
                row_obj['target_variant_len'] = 1
                row_obj['mutation_type'] = 'snp'
                row_obj['ref_state'] = ref_base
                row_obj['alt_state'] = alt_base
                state = 'ref'
                if base == alt_base and base != ref_base:
                    state = 'alt'
                row_obj['state'] = state
                row_obj['kseq'] = ''
                row_obj['klen'] = 0
                row_obj['homopolymer_len'] = 0
                row_obj['is_ambig_ok'] = True
                row_obj['is_kmer_found'] = True
                row_obj['is_kmer_length_ok'] = True
                row_obj['is_kmer_unique'] = True
                row_obj['is_valid'] = True

                if gene_feature is not None:
                    row_obj['gene_start'] = gene_start + 1
                    row_obj['gene_end'] = gene_end + 1
                    row_obj['cds_start'] = codon_var_start + 1
                    row_obj['cds_end'] = codon_var_end
                    row_obj['aa_name'] = aa_name
                    row_obj['aa_start'] = aa_var_start + 1
                    row_obj['aa_end'] = aa_var_start + 1

                row_obj['is_cds'] = is_cds
                row_obj['is_frame_shift'] = False
                row_obj['is_silent'] = is_silent
                row_obj['is_cds'] = is_cds
                row_obj['is_frame_shift'] = False
                row_obj['is_silent'] = is_silent
                row_obj['positive_genotypes'] = ",".join(sorted(rules[pos][base]['positive_genotypes']))
                row_obj['partial_genotypes'] = ",".join(sorted(rules[pos][base]['partial_genotypes']))

                counts = [0] * num_genotypes
                geno_found = set(rules[pos][base]['partial_genotypes']) | set(rules[pos][base]['positive_genotypes'])
                num_found = len(geno_found)

                for j in range(0, num_found):
                    counts[j] = 1
                row_obj['kmer_entropy'] = calc_shanon_entropy(counts)
                if max_entropy < row_obj['kmer_entropy']:
                    max_entropy = row_obj['kmer_entropy']
                scheme.append(row_obj)
                row_index+=1

    return scheme

def format_biohansel_scheme(biohansel_kmers,clade_obj):
    """
    A method to format the passed kmers and clade object to create a BioHansel
    compatible scheme elements.  Please refer to the file examples/small_test/cladeomatic/cladeomatic-biohansel.fasta
    for the scheme fasta file and the file examples/small_test/cladeomatic/cladeomatic-biohansel.meta.txt
    for the scheme metadata file.

    Parameters
    ----------
    biohansel_kmers : list
        The pre-formatted kmer list for use in the BioHansel scheme
    clade_obj : clade_worker object
        The clade_worker object (:class:`cladeomatic.clades.clade_worker`) that contains the clade data, selected genotypes, scheme data, and rule set for scheme creation

    Returns
    -------
    dict
        A dictionary containing the elements required for writing of a BioHansel formatted scheme - the fasta elements and the metadata elements.

    Notes
    -----
    Refer to https://github.com/phac-nml/biohansel for more thorough BioHansel documentation

    """
    valid_nodes = clade_obj.get_valid_nodes()
    node_heirarchy = {}
    delim = clade_obj.delim
    for sample_id in clade_obj.selected_genotypes:
        genotype = clade_obj.selected_genotypes[sample_id].split(delim)
        for i in range(0,len(genotype)):
            node_id = genotype[i]
            h = ".".join([str(x) for x in genotype[0:i+1]])
            node_heirarchy[node_id] = h

    clade_data = clade_obj.clade_data
    scheme = {'fasta':{},'meta':{}}
    for clade_id in valid_nodes:
        for i in range(0,len(clade_data[clade_id]['pos'])):
            p = clade_data[clade_id]['pos'][i]
            b = clade_data[clade_id]['bases'][i]
            if p in biohansel_kmers:
                start = biohansel_kmers[p][b]['start']
                if biohansel_kmers[p][b]['positive'] != '':
                    pid = "{}-{}".format(start,node_heirarchy[clade_id])
                    scheme['meta'][pid] = {'pos':p,'target_base':b}
                    scheme['fasta'][pid] = biohansel_kmers[p][b]['positive']
                if biohansel_kmers[p][b]['negative'] != '':
                    nid = "negative{}-{}".format(start, node_heirarchy[clade_id])
                    scheme['fasta'][nid] = biohansel_kmers[p][b]['negative']
                    scheme['meta'][nid] = {'pos': p, 'target_base': b}

    return scheme

def write_biohansel_scheme(scheme,fasta_file):
    """
    A method to format the passed scheme for BioHansel scheme writing.  Please refer to the file examples/small_test/cladeomatic/cladeomatic-biohansel.fasta
    for the scheme fasta file.

    Parameters
    ----------
    scheme : dict
        The BioHansel scheme data for writing to the fast file
    fasta_file : str
        The filepath for writing the BioHansel fasta file

    Notes
    -----
    Refer to https://github.com/phac-nml/biohansel for more thorough BioHansel documentation
    """
    fh = open(fasta_file,'w')
    for id in scheme:
        fh.write(">{}\n{}\n".format(id,scheme[id]))
    fh.close()

def write_biohansel_meta(scheme,out_file):
    """
    Method to write the BioHansel metadata table that includes the split
    kmers and the actual SNPs interrogated.  Please refer to the file examples/small_test/cladeomatic/cladeomatic-biohansel.meta.txt
    for the scheme metadata file.

    Parameters
    ----------
    scheme : dict
        The dictionary of scheme data for the metadata file output
    out_file : str
        The filepath for writing the BioHansel metadata file

    Notes
    -----
    Refer to https://github.com/phac-nml/biohansel for more thorough BioHansel documentation

    """
    fh = open(out_file,'w')
    fh.write("kmername\ttarget_pos\ttarget_base\n")
    for id in scheme:
        fh.write("{}\t{}\t{}\n".format(id,scheme[id]['pos'],scheme[id]['target_base']))
    fh.close()

def create_alt_psedo_sequence(ref_seq,positions,msa_base_counts,outfile):
    """
    This method creates a pseudo alternate sequence for the reference sequence
    with the SNP substitutions implemented.  A resulting fasta file is written for
    downstream use, refer to the file the file examples/small_test/cladeomatic/cladeomatic-altseq.fasta
    for more details.

    Parameters
    ----------
    ref_seq : str
        The reference sequence for processing
    positions : list
        A list of integer positions for the variants (SNPs)
    outfile : str
        The file path for writing the alternate sequence

    """

    #cast the string to an iterable list
    alt_seq = list(ref_seq)
    bases = ['A','T','C','G']
    basesr = range(0,len(bases))
    #loop through the positions of the mutations to sub in the alternate base
    for pos in positions:
        bases_present = set()
        for b in msa_base_counts[pos]:
            if msa_base_counts[pos][b] > 1:
                bases_present.add(b)
        candidates = set(bases) - bases_present
        if len(candidates) > 0:
            alt_seq[pos] = list(candidates)[0]

    alt_seq = ''.join(alt_seq)
    fh = open(outfile,'w')
    fh.write(">alt_seq\n{}\n".format(alt_seq))
    fh.close()

def run():
    """
    The main method to read the command line arguments and create the various output
    files as described in the readme.  Reads the tree or group file, reference file
    vcf file and metadata file.  Both the tree and group mode will validate the input
    files, label the nodes, create the distance matrix for the identified SNPs, filter
    for the relevant SNPs in the VCF, find the relevant kmers within the sequences
    submitted, determine the SNP scheme with and without kmers, and finally output
    all the files as per the documentation.  Please refer to :ref:`quickstart` for details on
    all output files.
    """

    cmd_args = parse_args()
    tree_file = cmd_args.in_nwk
    group_file = cmd_args.in_groups
    variant_file = cmd_args.in_var
    metadata_file = cmd_args.in_meta
    reference_file = cmd_args.reference
    prefix = cmd_args.prefix
    outdir = cmd_args.outdir
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
    max_site_ambig = cmd_args.max_site_ambig

    #initialize logging
    logging = init_console_logger(3)

    # Initialize Ray components
    os.environ['RAY_worker_register_timeout_seconds'] = '60'
    num_cpus = psutil.cpu_count(logical=False)
    if num_threads > num_cpus:
        num_threads = num_cpus

    if not ray.is_initialized():
        ray.init(num_cpus=num_threads)

    # Initialize output directory
    if not os.path.isdir(outdir):
        logging.info("Creating analysis directory {}".format(outdir))
        os.mkdir(outdir, 0o755)
    elif not force:
        logging.info("Error directory {} already exists, if you want to overwrite existing results then specify --force".format(outdir))
        sys.exit()

    #create the parameter log file of command line argument input for troubleshooting
    print_params(cmd_args, os.path.join(outdir,"{}-params.log".format(prefix)))

    #create the temporary processing directory
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

    # Validate input files for processing
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
    #parse the genbank or fasta reference sequences
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
    #parse the VCF file
    vcf_samples = set(vcfReader(variant_file).samples)
    logging.info("Reading metadata file")
    #parse the metadata file
    metadata = parse_metadata(metadata_file)
    #if the file is a tree file, validate and create a programmatic
    #ete3 tree file as per the rooting method passed - default is outgroup
    root_method = 'outgroup'
    if mode == 'tree':
        root_method = 'outgroup'

        # validate samples present in all input files
        #construct the tree using the midpoint rooting method
        ete_tree_obj = parse_tree(tree_file, logging, ete_format=1, set_root=True, resolve_polytomy=True,
                                  ladderize=True,
                                  method='midpoint')
        tree_samples = set(ete_tree_obj.get_leaf_names())

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

        sample_set = tree_samples | vcf_samples | set(metadata.keys())
        missing_samples = sample_set - tree_samples
    else:
        #parse the group file
        group_data = parse_group_file(group_file, delim=delim)
        group_samples = set(group_data['sample_map'].keys())
        sample_set = group_samples | vcf_samples | set(metadata.keys())
        missing_samples = sample_set - group_samples

    # Validate Sample ids across all inputs
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
        logging.error("Error specified minimum members to {} but there are only {} samples".format(min_member_count,
                                                                                                  len(sample_set)))
        sys.exit()

    #Parse all of the group information
    if mode == 'tree':
        group_data = parse_tree_groups(ete_tree_obj)
    else:
        group_data = parse_group_file(group_file, delim=delim)


    logging.info("Recreating fasta sequences from vcf")
    pseudo_seq_file = os.path.join(outdir, "pseudo.seqs.fasta")
    create_pseudoseqs_from_vcf(ref_seq_id,ref_seq[ref_seq_id], variant_file, pseudo_seq_file)

    #calculate distance matrix
    logging.info("Calculating SNP distance matrix")
    distance_matrix_file = os.path.join(outdir,"{}-dist.mat.txt".format(prefix))
    (stdout,stderr) = run_snpdists(pseudo_seq_file, distance_matrix_file, num_threads)
    if not os.path.isfile(distance_matrix_file) or os.path.getsize(distance_matrix_file) < MIN_FILE_SIZE:
        logging.error("snp-dists failed to produce a distance matrix, check the error message and try again:\n{}".format(stderr))
        sys.exit()

    #Filter vcf of invalid sites
    logging.info("Filtering VCF to remove sites which have > {} variants per site or > {} missing".format(max_states,max_site_ambig))
    filtered_vcf = os.path.join(outdir,"{}-filtered.vcf".format(prefix))
    count_snps_removed = filter_vcf(variant_file, filtered_vcf, max_states, max_site_ambig)
    logging.info("Removed {} SNPs from analysis and filtered vcf written to {}".format(count_snps_removed,filtered_vcf))

    #perform clade-snp work
    perform_compression = True
    if no_compression:
        perform_compression = False
    logging.info("Performing canonical SNP detection")
    cw = clade_worker(filtered_vcf, metadata , distance_matrix_file, group_data, ref_seq[ref_seq_id], mode,perform_compression=perform_compression,delim=delim,
                      min_snp_count=min_snp_count, max_snps=max_snp_count, max_states=max_states, min_members=min_member_count,
                 min_inter_clade_dist=1, num_threads=num_threads,rcor_thresh=rcor_thresh)
    snp_genotype_rules = cw.get_genotype_snp_rules()
    logging.info("Read {} variant positions from {}".format(cw.num_positions,filtered_vcf))
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
    kw = kmer_worker(ref_seq[ref_seq_id], pseudo_seq_file, analysis_dir, prefix, klen, genotype_map, snp_genotype_rules,max_ambig=max_ambig, min_perc=min_perc,
                 target_positions=target_positions, num_threads=num_threads)

    kmers = kw.extracted_kmers
    #write the extracted kmer file with kmers determined above
    fh = open(os.path.join(outdir,"{}-extracted.kmers.txt".format(prefix)),'w')
    fh.write("kseq\taln_start\taln_end\tis_valid\ttarget_positions\tgenotype\tgenotype_count\n")
    for kseq in kmers:
        for genotype in kmers[kseq]['genotype_counts']:
            count = kmers[kseq]['genotype_counts'][genotype]
            row = [
                kseq,
                kmers[kseq]['aln_start'],
                kmers[kseq]['aln_end'],
                kmers[kseq]['is_valid'],
                ";".join([str(x) for x in kmers[kseq]['target_positions']]),
                genotype,
                count
            ]
            fh.write("{}\n".format("\t".join([str(x) for x in row])))
    fh.close()


    #Update based on positions which could not be assigned a kmer
    positions_missing_kmer = kw.positions_missing_kmer
    if len(positions_missing_kmer) > 0:
        cw.remove_snps(list(positions_missing_kmer.keys()))
        cw.update()
    logging.info("A total of {} genotyping positions removed due to no valid kmers found: {}".format(len(positions_missing_kmer),positions_missing_kmer))
    logging.info("Final set of {} genotyping positions selected".format(len(cw.selected_positions)))
    write_genotypes(cw.selected_genotypes, os.path.join(outdir, "{}-genotypes.selected.txt".format(prefix)))

    logging.info("Writing distance based clustering")

    dist_header = "{}\n".format("\t".join(["sample_id"] + ["thresh-{}".format(";".join([str(x) for x in cw.dist_thresholds]))]))
    dist_nomenclature = {}
    for sample_id in cw.sample_linkage_clusters:
        genotype = cw.delim.join([str(x) for x in cw.sample_linkage_clusters[sample_id]])
        dist_nomenclature[sample_id] = genotype
    write_genotypes(dist_nomenclature, os.path.join(outdir, "{}-genotypes.distance.txt".format(prefix)),header=dist_header)

    #if sepcified, prune the scheme for the maximum number of snps per node
    if max_snp_count > 1:
        cw.prune_snps()
        variant_pos = set(cw.variant_positions)
        selected_positions = set(cw.get_selected_positions())
        kw.remove_scheme_pos(variant_pos - selected_positions)

    logging.info("Creating schemes")
    if len(ref_features) > 0:
        scheme = create_scheme(SCHEME_HEADER,ref_features,kw,cw.selected_genotypes,trans_table=11)
        snp_scheme = create_snp_scheme(SCHEME_HEADER,ref_features,cw,trans_table=11)
    else:
        scheme = create_scheme(SCHEME_HEADER, ref_seq, kw, cw.selected_genotypes, trans_table=11)
        snp_scheme = create_snp_scheme(SCHEME_HEADER, ref_seq, cw, trans_table=11)

    write_scheme(SCHEME_HEADER, scheme, os.path.join(outdir, "{}-kmer.scheme.txt".format(prefix)))
    write_scheme(SCHEME_HEADER, snp_scheme, os.path.join(outdir, "{}-snp.scheme.txt".format(prefix)))

    #creating the BioHansel files
    bh_data = format_biohansel_scheme(kw.biohansel_kmers, cw)
    write_biohansel_scheme(bh_data['fasta'], os.path.join(outdir, "{}-biohansel.fasta".format(prefix)))
    write_biohansel_meta(bh_data['meta'], os.path.join(outdir, "{}-biohansel.meta.txt".format(prefix)))

    if not keep_tmp:
        logging.info("Removing temporary analysis folder: {}".format(analysis_dir))
        shutil.rmtree(analysis_dir)

    msa_base_counts = kw.msa_base_counts
    create_alt_psedo_sequence(ref_seq[ref_seq_id], cw.selected_positions, msa_base_counts, os.path.join(outdir,"{}-altseq.fasta".format(prefix)))

    logging.info("Analysis complete")
