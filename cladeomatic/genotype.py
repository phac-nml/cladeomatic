import sys
import os
import ray
from cladeomatic.version import __version__
from cladeomatic.utils.vcfhelper import vcfReader
from cladeomatic.utils import init_console_logger
from cladeomatic.utils import parse_metadata
import pandas as pd
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from cladeomatic.constants import GENOTYPE_REPORT_HEADER, MIN_FILE_SIZE
from datetime import datetime


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
    parser.add_argument('--sample_meta', type=str, required=True, help='Tab delimited sample metadata', default=None)
    parser.add_argument('--genotype_meta', type=str, required=False, help='Tab delimited genotype metadata', default=None)
    parser.add_argument('--outfile', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('--max_missing_positions', type=int, required=False, help='Maximum number of missing positions for the genotype', default=1)
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

def parse_scheme_features(scheme_file):

    features = {
        'scheme':{},
        'num_positions':0,
        'num_mutations':0,
        'total_scheme_features':0,
        'mutation_lookup':{}
    }

    scheme = {}
    df = pd.read_csv(scheme_file, sep="\t", header=0, low_memory=False)
    snp_states = {}
    num_positions = df['variant_start'].nunique()
    num_features = df['mutation_key'].nunique()
    features['total_scheme_features'] = len(df)
    features['num_positions'] = num_positions
    features['num_mutations'] = num_features

    for row in df.itertuples():
        state = row.state
        mutation_key = row.mutation_key
        target_variant = row.target_variant
        variant_start = int(row.variant_start)
        if not variant_start in features['mutation_lookup']:
            features['mutation_lookup'][variant_start] = {'alt':{},'ref':{}}
        features['mutation_lookup'][variant_start][state][target_variant] = mutation_key

        positive_genotypes = row.positive_genotypes
        if isinstance(positive_genotypes, float):
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')

        genotypes = positive_genotypes
        for genotype in positive_genotypes:

            if not genotype in scheme:
                scheme[genotype] = {'positive': {}, 'partial': {}, 'allowed': {}}
            if not variant_start in scheme[genotype]['positive']:
                scheme[genotype]['positive'][variant_start] = set()
            scheme[genotype]['positive'][variant_start].add(target_variant)

        partial_genotypes = row.partial_genotypes

        if not variant_start in snp_states:
            snp_states[variant_start] = set()
        snp_states[variant_start].add(target_variant)

        if isinstance(partial_genotypes, float):
            partial_genotypes = []
        else:
            partial_genotypes = partial_genotypes.split(',')
        genotypes += partial_genotypes
        if len(partial_genotypes) == 0:
            continue
        for genotype in partial_genotypes:
            if not genotype in scheme:
                scheme[genotype] = {'positive': {}, 'partial': {}, 'allowed': {}}
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

    features['scheme'] = scheme
    return features



def parse_scheme_genotypes(scheme_file):
    scheme = {}
    df = pd.read_csv(scheme_file, sep="\t", header=0, low_memory=False)
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

    return scheme

@ray.remote
def call_genotypes(genotype_rules,metadata,variants,max_dist=0):
    result = {}
    for sample_id in metadata:

        if not sample_id in variants:
            continue

        result[sample_id] = {
            'predicted_genotype(s)':[],
            'predicted_genotype_dist': [],
            'genoytpe_results':[],
            'genoytpe_dists': {},

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
                result[sample_id]['predicted_genotype_dist'].append(dist)
                result[sample_id]['genoytpe_results'].append({'match':genoytpe_results[genotype]['match'],
                                    'mismatch':genoytpe_results[genotype]['mismatch']})
                pdist = dist

        del result[sample_id]['genoytpe_dists']

    return result


def convert_features_to_mutations(feature_lookup,snp_profile):
    found_features = {
        'alt':set(),
        'ref':set()
    }
    for pos in snp_profile:
        if pos not in feature_lookup:
            continue
        base = snp_profile[pos]
        for state in feature_lookup[pos]:
            if base in feature_lookup[pos][state]:
                found_features[state].add(feature_lookup[pos][state][base])
                break


    return found_features





def write_genotype_calls(header,scheme_name,outfile,genotype_results,sample_metadata,genotype_meta, scheme_data, sample_variants,min_positions=1):

    #get all additional_fields
    sample_fields = set()
    for sample_id in sample_metadata:
        for field in sample_metadata[sample_id]:
            sample_fields.add(field)
    sample_fields = sorted(list(sample_fields))

    genotype_fields = set()
    for genotype in genotype_meta:
        for field in genotype_meta[genotype]:
            genotype_fields.add(field)
    genotype_fields = sorted(list(genotype_fields))
    header = header + sample_fields + genotype_fields

    #initialize file
    fh = open(outfile, 'w')
    fh.write("{}\n".format("\t".join([str(x) for x in header])))

    analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")


    num_positions = scheme_data['num_positions']
    total_scheme_features = scheme_data['total_scheme_features']


    for sample_id in genotype_results:
        row = {}
        status = 'Pass'
        for field_id in header:
            row[field_id] = ''
        row['sample_id'] = sample_id
        row['predicted_genotype'] = '-'
        row['scheme'] = scheme_name
        row['analysis_date'] = analysis_date
        row['total_scheme_features'] = total_scheme_features
        row['unique_scheme_positions'] = num_positions

        detected_mutations = convert_features_to_mutations(scheme_data['mutation_lookup'],sample_variants[sample_id])
        row['num_detected_positions'] =  len(detected_mutations['alt'] | detected_mutations['ref'])
        row['detected_alt_mutations'] = ";".join([str(x) for x in sorted(list(detected_mutations['alt']))])

        for field_id in sample_metadata[sample_id]:
            row[field_id] = sample_metadata[sample_id][field_id]

        if len(genotype_results[sample_id]['predicted_genotype(s)']) == 1:
            row['predicted_genotype'] = genotype_results[sample_id]['predicted_genotype(s)'][0]
            row['predicted_genotype_distance'] = genotype_results[sample_id]['predicted_genotype_dist'][0]
            genotype = row['predicted_genotype']
            if genotype in genotype_meta:
                for field_id in genotype_meta[genotype]:
                    row[field_id] = genotype_meta[genotype][field_id]
        elif len(genotype_results[sample_id]['predicted_genotype(s)']) > 1:
            status = 'Warning'
            row['qc_messages'] = "Ambiguous genotype assignement, possible genotypes: {}".format(";".join([str(x) for x in genotype_results[sample_id]['predicted_genotype(s)']]))
        else:
            status = 'Warning'
            row['qc_messages'] = 'No genotypes were compatible with the sample'


        if row['num_detected_positions'] < min_positions:
            status = 'Fail'
            row['qc_messages'] = 'Sample is missing too many positions for genotype call'

        row['qc_status'] = status
        fh.write("{}\n".format("\t".join([str(x) for x in row.values()])))

    fh.close()

def is_valid_file(filename):
    status = True
    if not os.path.isfile(filename) or os.path.getsize(filename) <= MIN_FILE_SIZE :
        status = False
    return status


def run():
    cmd_args = parse_args()
    scheme_file = cmd_args.in_scheme
    variant_file = cmd_args.in_var
    metadata_file = cmd_args.sample_meta
    genotype_meta_file = cmd_args.genotype_meta
    outfile = cmd_args.outfile
    num_threads = cmd_args.num_threads
    max_missing_positions = cmd_args.max_missing_positions

    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True, num_cpus=num_threads)


    logging = init_console_logger(3)
    logging.info("Starting analysis")

    logging.info("Reading metadata file {}".format(metadata_file))
    if not is_valid_file(metadata_file):
        logging.error("Error file {} was not found or is empty".format(metadata_file))
        sys.exit()
    sample_metadata = parse_metadata(metadata_file)
    sample_meta_fields = list(sample_metadata[list(sample_metadata.keys())[0]].keys())

    genotype_metadata = {}
    if genotype_meta_file is not None:
        if not is_valid_file(genotype_meta_file):
            logging.error("Error file {} was not found or is empty".format(genotype_meta_file))
            sys.exit()
        logging.info("Reading metadata file {}".format(genotype_meta_file))
        genotype_metadata = parse_metadata(genotype_meta_file)


    logging.info("Reading scheme file {}".format(scheme_file))
    if not is_valid_file(scheme_file):
        logging.error("Error file {} was not found or is empty".format(scheme_file))
        sys.exit()

    scheme_data = parse_scheme_features(scheme_file)
    genotype_rules = scheme_data['scheme']

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
    min_positions = len(valid_positions) -max_missing_positions

    logging.info("Extracted {} genotyping positions".format(len(valid_positions)))
    logging.info("Reading snp data from vcf file {}".format(variant_file))
    variants = get_snp_profiles(valid_positions, variant_file)


    logging.info("Calling genotypes for {} samples based on {} genotypes".format(len(sample_metadata),num_genotypes))
    batch_size = int(len(sample_metadata) / num_threads)
    ray_results = []
    if batch_size < 1:
        batch_size = 1
    sub_metadata = {}
    sub_variants = {}
    for sample_id in sample_metadata:
        sub_metadata[sample_id] = sample_metadata[sample_id]
        sub_variants[sample_id] = variants[sample_id]
        if len(sub_metadata) == batch_size:
            ray_results.append(call_genotypes.remote(rule_id, sub_metadata, sub_variants))
            sub_metadata = {}
            sub_variants = {}
    if len(sub_metadata) > 0:
        ray_results.append(call_genotypes.remote(rule_id, sub_metadata, sub_variants))
    results = ray.get(ray_results)

    genotype_results = {}
    for r in results:
        for sample_id in r:
            genotype_results[sample_id] = r[sample_id]
    del(results)
    del(ray_results)

    write_genotype_calls(GENOTYPE_REPORT_HEADER, os.path.basename(scheme_file), outfile, genotype_results, sample_metadata, genotype_metadata, scheme_data,variants, min_positions)

    logging.info("Analysis complete")

