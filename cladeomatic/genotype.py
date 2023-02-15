import sys

import ray
from cladeomatic.version import __version__
from cladeomatic.utils.vcfhelper import vcfReader
from cladeomatic.utils import init_console_logger
from cladeomatic.utils import parse_metadata
import pandas as pd
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)


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
    parser.add_argument('--in_meta', type=str, required=True, help='Tab delimited genotype metadata', default=None)
    parser.add_argument('--outfile', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('--num_threads', type=int, required=False, help='Number of threads to use', default=1)
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
    snp_states = {}
    for row in df.itertuples():
        target_variant = row.target_variant
        variant_start = int(row.variant_start)
        positive_genotypes = row.positive_genotypes
        if isinstance(positive_genotypes,float):
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')

        genotypes = positive_genotypes
        for genotype in positive_genotypes:

            if not genotype in scheme:
                scheme[genotype] = {'positive':{},'partial':{},'allowed':{}}
            if not variant_start in scheme[genotype]['positive']:
                scheme[genotype]['positive'][variant_start] = set()
            scheme[genotype]['positive'][variant_start].add(target_variant)

        partial_genotypes = row.partial_genotypes

        if not variant_start in snp_states:
            snp_states[variant_start] = set()
        snp_states[variant_start].add(target_variant)

        if isinstance(partial_genotypes,float):
            partial_genotypes = []
        else:
            partial_genotypes = partial_genotypes.split(',')
        genotypes+= partial_genotypes
        if len(partial_genotypes) == 0:
            continue
        for genotype in partial_genotypes:
            if not genotype in scheme:
                scheme[genotype] = {'positive':{},'partial':{},'allowed':{}}
            if not variant_start in scheme[genotype]['partial']:
                scheme[genotype]['partial'][variant_start] = set()
            scheme[genotype]['partial'][variant_start].add(target_variant)


    for genotype in scheme:
        for pos in snp_states:
            allowed_bases = set()
            if pos in scheme[genotype]['partial']:
                allowed_bases = allowed_bases | set(scheme[genotype]['partial'][pos])
            if pos in scheme[genotype]['positive']:
                allowed_bases = allowed_bases | set(scheme[genotype]['positive'][pos])
            if len(allowed_bases) == 0:
                allowed_bases = snp_states[pos]
            scheme[genotype]['allowed'][pos] = allowed_bases


    genotypes = list(scheme.keys())

    '''
        for i in range(0,len(genotypes)):
        g1 = genotypes[i]
        a1 = scheme[g1]['allowed']
        for k in range(i+1,len(genotypes)):
            g2 = genotypes[k]
            dist = 0
            a2 = scheme[g2]['allowed']
            for pos in a1:
                int1 = a1[pos] & a2[pos]
                if len(int1) == 1:
                    dist+=1
    '''



    return scheme

@ray.remote
def call_genotypes(genotype_rules,metadata,variants,max_dist=0):
    result = {}
    for sample_id in metadata:
        if not 'genotype' in metadata[sample_id]:
            continue
        if not sample_id in variants:
            continue

        result[sample_id] = {
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
            found_base = set(variants[sample_id][pos])
            if '*' in found_base  or '-' in found_base or 'N' in found_base:
                continue

            for genotype in genotype_rules:
                allowed_bases = genotype_rules[genotype]['allowed'][pos]

                if len(found_base & allowed_bases) == 1:
                    genoytpe_results[genotype]['match'][pos] = found_base

                else:
                    genoytpe_results[genotype]['mismatch'][pos] = found_base

        for genotype in genoytpe_results:
            matches = len(genoytpe_results[genotype]['match'])
            mismatches = len(genoytpe_results[genotype]['mismatch'])
            total = matches + mismatches
            if total > 0:
                dists[genotype] = 1 - matches /total

                
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
        del result[sample_id]['genoytpe_dists']

    return result

def run():
    cmd_args = parse_args()
    scheme_file = cmd_args.in_scheme
    variant_file = cmd_args.in_var
    metadata_file = cmd_args.in_meta
    outfile = cmd_args.outfile
    num_threads = cmd_args.num_threads

    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True, num_cpus=num_threads)


    logging = init_console_logger(3)
    logging.info("Starting analysis")

    logging.info("Reading metadata file {}".format(metadata_file))
    metadata = parse_metadata(metadata_file)
    fields = list(metadata[list(metadata.keys())[0]].keys())



    logging.info("Reading scheme file {}".format(scheme_file))
    genotype_rules = parse_scheme_genotypes(scheme_file)
    for genotype in genotype_rules:
        if len(genotype_rules[genotype]['positive']) == 0:
            logging.warn("Genotype {} has no required kmers".format(genotype))



    rule_id = ray.put(genotype_rules)
    num_genotypes = len(genotype_rules)


    valid_positions = []
    for genotype in genotype_rules:
        for state in genotype_rules[genotype]:
            valid_positions += list(genotype_rules[genotype][state].keys())
    valid_positions = list(set(valid_positions))

    logging.info("Extracted {} genotyping positions".format(len(valid_positions)))

    logging.info("Reading snp data from vcf file {}".format(variant_file))
    variants = get_snp_profiles(valid_positions, variant_file)

    logging.info("Calling genotypes for {} samples based on {} genotypes".format(len(metadata),num_genotypes))
    batch_size = int(len(metadata) / num_threads)
    ray_results = []
    if batch_size < 1:
        batch_size = 1
    sub_metadata = {}
    sub_variants = {}
    for sample_id in metadata:
        sub_metadata[sample_id] = metadata[sample_id]
        sub_variants[sample_id] = variants[sample_id]
        if len(sub_metadata) == batch_size:
            ray_results.append(call_genotypes.remote(rule_id, sub_metadata, sub_variants))
            sub_metadata = {}
            sub_variants = {}
    if len(sub_metadata) > 0:
        ray_results.append(call_genotypes.remote(rule_id, sub_metadata, sub_variants))
    results = ray.get(ray_results)
    genoytpe_results = {}
    for r in results:
        for sample_id in r:
            genoytpe_results[sample_id] = r[sample_id]
    del(results)
    del(ray_results)

    fh = open(outfile,'w')
    header = ['sampleid','predicted_genotype']
    for field in fields:
        header.append(field)

    fh.write("{}\n".format("\t".join([str(x) for x in header])))
    num_fields = len(fields)
    for sample_id in genoytpe_results:
        genotypes = [str(x) for x in genoytpe_results[sample_id]['predicted_genotype(s)']]
        row = [
            sample_id,
            ",".join(genotypes),
        ]
        for i in range(0,num_fields):
            field_name = fields[i]
            contents = []
            for g in genotypes:
                value = ''
                if g in metadata and field_name in metadata[g]:
                    value = metadata[g][field_name]
                contents.append(value)
            row.append(";".join([str(x) for x in contents]))
        fh.write("{}\n".format("\t".join(row)))

    fh.close()

    logging.info("Analysis complete")

