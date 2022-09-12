import os
import sys

from cladeomatic.version import __version__
from cladeomatic.utils.vcfhelper import vcfReader
from cladeomatic.utils import init_console_logger
from cladeomatic.utils import parse_metadata
import pandas as pd
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from sklearn.metrics import f1_score

def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Clade-O-Matic: Genotyping scheme development v. {}".format(__version__),
        formatter_class=CustomFormatter)
    parser.add_argument('--in_var', type=str, required=True,
                        help='Either Variant Call SNP data (.vcf) or TSV SNP data (.txt)')
    parser.add_argument('--in_scheme', type=str, required=True, help='Tab delimited scheme file produced by clade-o-matic',
                        default=None)
    parser.add_argument('--in_meta', type=str, required=True, help='Tab delimited file of genotype assignments', default=None)
    parser.add_argument('--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('--prefix', type=str, required=False, help='Prefix for output files', default='cladeomatic')
    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()


def get_snp_profiles(valid_positions, vcf_file):
    '''
    Accepts SNP position list and vcf file
    :param valid_positions: list of integers
    :param vcf_file: str path to vcf or tsv snp data
    :return: dict of snp_data data structure
    '''
    vcf = vcfReader(vcf_file)
    data = vcf.process_row()

    samples = vcf.samples
    profiles = {}
    count_snps = 0
    if data is not None:
        count_snps += 1
    for sample_id in samples:
        profiles[sample_id] = {}
    if data is None:
        return profiles
    while data is not None:
        pos = int(data['POS'])
        if not pos in valid_positions:
            data = vcf.process_row()
            continue
        for sample_id in samples:
            base = data[sample_id]
            profiles[sample_id][pos] = base
        count_snps += 1

        data = vcf.process_row()

    return profiles

def parse_scheme_genotypes(scheme_file):
    scheme = {}
    df = pd.read_csv(scheme_file, sep="\t", header=0, low_memory=False)
    variant_positions = list(df['variant_start'].unique())
    geno_seqs = {}
    for row in df.itertuples():
        target_variant = row.target_variant
        variant_start = int(row.variant_start)
        positive_genotypes = row.positive_genotypes
        if isinstance(positive_genotypes,float):
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')
        if len(positive_genotypes) == 0:
            continue
        genotypes = positive_genotypes
        for genotype in positive_genotypes:

            if not genotype in scheme:
                scheme[genotype] = {'positive':{},'partial':{}}
            scheme[genotype]['positive'][variant_start] = target_variant

        partial_genotypes = row.partial_genotypes

        if isinstance(partial_genotypes,float):
            partial_genotypes = []
        else:
            partial_genotypes = partial_genotypes.split(',')
        genotypes+= partial_genotypes
        if len(partial_genotypes) == 0:
            continue
        for genotype in partial_genotypes:
            if not genotype in scheme:
                scheme[genotype] = {'positive':{},'partial':{}}
            scheme[genotype]['partial'][variant_start] = target_variant

        for genotype in genotypes:
            if not genotype in geno_seqs:
                geno_seqs[genotype] = {}
                for pos in variant_positions:
                    geno_seqs[genotype][pos] = 'N'
            if geno_seqs[genotype][variant_start] == 'N':
                geno_seqs[genotype][variant_start] = target_variant
            else:
                if geno_seqs[genotype][variant_start] != target_variant:
                    geno_seqs[genotype][variant_start] = "N"
    return scheme

def call_genotypes(genotype_rules,metadata,variants,max_dist=0,n_threads=1):
    result = {}
    for sample_id in metadata:
        if not 'genotype' in metadata[sample_id]:
            continue
        if not sample_id in variants:
            continue
        genotype = metadata[sample_id]['genotype']
        result[sample_id] = {
            'submitted_genotype':genotype,
            'predicted_genotype(s)':[],
            'predicted_genotype_dist': 1,
            'genoytpe_results':{},
            'genoytpe_dists': {}
        }
        genoytpe_results = {}
        dists = {}
        for genotype in genotype_rules:
            genoytpe_results[genotype] = {'match':{},'mismatch':{}}
            dists[genotype] = 1

        for pos in variants[sample_id]:
            found_base = variants[sample_id][pos]
            for genotype in genotype_rules:
                if pos in genotype_rules[genotype]['positive']:
                    target_base = genotype_rules[genotype]['positive'][pos]
                    if found_base == target_base:
                        genoytpe_results[genotype]['match'][pos] = found_base
                    else:
                        genoytpe_results[genotype]['mismatch'][pos] = found_base
                elif pos in genotype_rules[genotype]['partial']:
                    target_base = genotype_rules[genotype]['partial'][pos]
                    if found_base == target_base:
                        genoytpe_results[genotype]['match'][pos] = found_base

        for genotype in genoytpe_results:
            matches = len(genoytpe_results[genotype]['match'])
            mismatches = len(genoytpe_results[genotype]['mismatch'])
            total = matches + mismatches
            if total > 0:
                dists[genotype] = 1 - matches /total

        result[sample_id]['genoytpe_results'] = genoytpe_results
        result[sample_id]['genoytpe_dists'] =  {k: v for k, v in sorted(dists.items(), key=lambda item: item[1])}
        pdist = 1
        for genotype in result[sample_id]['genoytpe_dists']:
            dist =  result[sample_id]['genoytpe_dists'][genotype]
            if dist <= pdist and dist <= max_dist:
                result[sample_id]['predicted_genotype(s)'].append(genotype)
                result[sample_id]['predicted_genotype_dist'] = dist
                pdist = dist

        num_geno = len(result[sample_id]['predicted_genotype(s)'])
        if num_geno > 1:
            filt = []
            for i in range(0,num_geno):
                is_substring = False
                for k in range(i+1,num_geno):
                    if "{}.".format(result[sample_id]['predicted_genotype(s)'][i]) in result[sample_id]['predicted_genotype(s)'][k]:
                        is_substring = True
                if not is_substring:
                    filt.append(result[sample_id]['predicted_genotype(s)'][i])
            result[sample_id]['predicted_genotype(s)'] = filt



    return result

def run():
    cmd_args = parse_args()
    scheme_file = cmd_args.in_scheme
    variant_file = cmd_args.in_var
    metadata_file = cmd_args.in_meta
    prefix = cmd_args.prefix
    outdir = cmd_args.outdir

    logging = init_console_logger(3)
    logging.info("Starting analysis")

    if not os.path.isdir(outdir):
        logging.info("Creating temporary analysis directory {}".format(outdir))
        os.mkdir(outdir, 0o755)

    logging.info("Reading metadata file {}".format(metadata_file))
    metadata = parse_metadata(metadata_file)

    logging.info("Reading scheme file {}".format(scheme_file))
    genotype_rules = parse_scheme_genotypes(scheme_file)


    valid_positions = []
    for genotype in genotype_rules:
        for state in genotype_rules[genotype]:
            valid_positions += list(genotype_rules[genotype][state].keys())
    valid_positions = list(set(valid_positions))

    logging.info("Extracted {} genotyping positions".format(len(valid_positions)))

    logging.info("Reading snp data from vcf file {}".format(variant_file))
    variants = get_snp_profiles(valid_positions, variant_file)

    logging.info("Calling genotypes for {} samples".format(len(metadata)))
    genoytpe_results = call_genotypes(genotype_rules, metadata, variants)

    fh = open(os.path.join(outdir,"{}-scheme.calls.txt".format(prefix)),'w')
    fh.write("sample_id\tsubmited_genotype\tpredicted_genotype\tis_match\n")
    for sample_id in genoytpe_results:
        g = genoytpe_results[sample_id]['submitted_genotype']
        c = ",".join([str(x) for x in genoytpe_results[sample_id]['predicted_genotype(s)']])
        m = g == c
        fh.write("{}\t{}\t{}\t{}\n".format(sample_id, g, c,m))

    logging.info("Calcualting F1 for {} samples".format(len(metadata)))
    truth = []
    pred = []
    group_samples = {}
    num_ambig = 0
    num_correct = 0
    for sample_id in genoytpe_results:
        genotype = str(genoytpe_results[sample_id]['submitted_genotype'])

        if genotype not in group_samples:
            group_samples[genotype] = {'truth': [], 'pred': []}
        pgenotypes = genoytpe_results[sample_id]['predicted_genotype(s)']

        if len(pgenotypes) > 1:
            num_ambig+=1

        if len(pgenotypes) >= 1 and genotype == str(pgenotypes[0]):
            pred.append(genotype)
            group_samples[genotype]['pred'].append(sample_id)
            num_correct+=1
        else:
            pred.append('X')
        truth.append(genotype)
        group_samples[genotype]['truth'].append(sample_id)

    scheme_f1 = f1_score(truth,pred,average='micro')
    fh = open(os.path.join(outdir,"{}-scheme.scores.txt".format(prefix)),'w')
    fh.write("genotype\tnum_true\tnum_pred\tf1\n")
    fh.write("overall\t{}\t{}\t{}\n".format(len(metadata),num_correct,scheme_f1))
    for genotype in group_samples:
        num_true = len(group_samples[genotype]['truth'])
        samples = list(set(group_samples[genotype]['truth'] + group_samples[genotype]['pred']))
        num_samples = len(samples)
        truth = [0] * num_samples
        pred = [0] * num_samples
        for i in range(0,num_samples):
            sample_id = samples[i]
            if sample_id in group_samples[genotype]['truth']:
                truth[i] = 1
            if sample_id in group_samples[genotype]['pred']:
                pred[i] = 1
        geno_f1 = f1_score(truth, pred)
        fh.write("{}\t{}\t{}\t{}\n".format(genotype,num_true,num_samples,geno_f1))
    fh.close()

    logging.info("Analysis complete")

