import os
import ray
from cladeomatic.version import __version__
from cladeomatic.utils import init_console_logger
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from sklearn.metrics import f1_score
from cladeomatic.genotype import parse_scheme_features, get_snp_profiles
from cladeomatic.constants import GENOTYPE_REPORT_HEADER, SCHEME_SCORES_REPORT_HEADER, SCHEME_SAMPLE_REPORT_HEADER
from cladeomatic.utils import parse_metadata
from datetime import datetime

def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Clade-O-Matic: Benchmarking Genotyping scheme development v. {}".format(__version__),
        formatter_class=CustomFormatter)
    parser.add_argument('--in_genotype', type=str, required=True,
                        help='Genotype report made by genotyper')
    parser.add_argument('--in_var', type=str, required=True,
                        help='Either Variant Call SNP data (.vcf) or TSV SNP data (.txt)')
    parser.add_argument('--in_scheme', type=str, required=True,
                        help='Tab delimited scheme file produced by clade-o-matic',
                        default=None)
    parser.add_argument('--submitted_genotype_col', type=str, required=True, help='Name of column containing submitted genotype')
    parser.add_argument('--predicted_genotype_col', type=str, required=True, help='Name of column containing predicted genotype')
    parser.add_argument('--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('--prefix', type=str, required=False, help='Output Directory to put results',default='cladeomatic')
    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()

def parse_genotype_report(header,file):
    return parse_metadata(file,column=header[0])

def convert_genotype_report(data_dict,predicted_field_name,submitted_field_name):
    p = []
    s = []
    samples = []
    for sample_id in data_dict:
        samples.append(str(sample_id))
        p.append(str(data_dict[sample_id][predicted_field_name]))
        s.append(str(data_dict[sample_id][submitted_field_name]))
    return {'sample_ids':samples,'predicted':p,'submitted':s}

def filter_genotypes_include(labels,name):
    p = []
    s = []
    samples = []
    for i in range(0,len(labels['sample_ids'])):
        if name in labels['predicted'][i] or name in labels['submitted'][i]:
            samples.append(labels['sample_ids'][i])
            p.append(labels['predicted'][i])
            s.append(labels['submitted'][i])
    return {'sample_ids':samples,'predicted':p,'submitted':s}

def filter_genotypes_exclude(labels,name):
    p = []
    s = []
    samples = []
    for i in range(0,len(labels['sample_ids'])):
        if name in labels['predicted'][i] or name in labels['submitted'][i]:
            continue
        samples.append(labels['sample_ids'][i])
        p.append(labels['predicted'][i])
        s.append(labels['submitted'][i])
    return {'sample_ids':samples,'predicted':p,'submitted':s}

def write_scheme_scores(header,file,scheme_name,scores):
    fh = open(file,'w')
    fh.write("{}\n".format("\t".join(header)))
    analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    for label in scores:
        row = {}
        for field in header:
            row[field] = ''
        row['label'] = label
        row['scheme'] = scheme_name
        row['analysis_date'] = analysis_date
        for field in scores[label]:
            row[field] = scores[label][field]
        fh.write("{}\n".format("\t".join([str(x) for x in row.values()])))
    fh.close()

def get_problem_bases(profile,rule):
    mismatches = {}
    for pos in profile:
        found_base = set(profile[pos])
        if '*' in found_base or '-' in found_base or 'N' in found_base:
            continue
        allowed_bases = rule['allowed'][pos]
        if len(found_base & allowed_bases) != 1:
            mismatches[pos] = found_base
    return mismatches

def write_updated_genotype_report(header,file,scheme_name,variants,data_dict,genotype_rules,predicted_field_name,submitted_field_name):
    fh = open(file,'w')
    fh.write("{}\n".format("\t".join(header)))
    analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    for sample_id in data_dict:
        p = data_dict[sample_id][predicted_field_name]
        s = data_dict[sample_id][submitted_field_name]
        is_match = p == s
        mismatches = {}
        if not is_match:
            mismatches = get_problem_bases(variants[sample_id],genotype_rules[s])
        row = {}
        for field in header:
            row[field] = ''
        row['sample_id'] = sample_id
        row['scheme'] = scheme_name
        row['analysis_date'] = analysis_date
        row['submitted_genotype'] = s
        row['predicted_genotype'] = p
        row['is_match'] = is_match
        m = []
        for pos in mismatches:
            m.append("{}:{}".format(pos,mismatches[pos]))
        row['problem_feature(s)'] = ",".join(m)
        fh.write("{}\n".format("\t".join([str(x) for x in row.values()])))
    fh.close()
    return

def get_genotype_counts(data_dict,field_name='genotype'):
    counts = {}
    for sample_id in data_dict:
        genotype = data_dict[sample_id][field_name]
        if not genotype in counts:
            counts[genotype] = 0
        counts[genotype]+=1
    return counts



def run():
    cmd_args = parse_args()
    scheme_file = cmd_args.in_scheme
    variant_file = cmd_args.in_var
    geno_file = cmd_args.in_genotype
    prefix = cmd_args.prefix
    outdir = cmd_args.outdir
    predicted_field_name = cmd_args.predicted_genotype_col
    submitted_field_name = cmd_args.submitted_genotype_col

    logging = init_console_logger(3)
    logging.info("Starting analysis")

    if not os.path.isdir(outdir):
        logging.info("Creating temporary analysis directory {}".format(outdir))
        os.mkdir(outdir, 0o755)

    logging.info("Reading scheme file {}".format(scheme_file))
    scheme_data = parse_scheme_features(scheme_file)
    genotype_rules = scheme_data['scheme']

    for genotype in genotype_rules:
        if len(genotype_rules[genotype]['positive']) == 0:
            logging.warn("Genotype {} has no required kmers".format(genotype))

    valid_positions = []
    for genotype in genotype_rules:
        for state in genotype_rules[genotype]:
            valid_positions += list(genotype_rules[genotype][state].keys())
    valid_positions = list(set(valid_positions))

    logging.info("Extracted {} genotyping positions".format(len(valid_positions)))

    logging.info("Reading snp data from vcf file {}".format(variant_file))
    variants = get_snp_profiles(valid_positions, variant_file)


    data_dict = parse_genotype_report(GENOTYPE_REPORT_HEADER,geno_file)
    num_samples = len(data_dict)
    labels = convert_genotype_report(data_dict, predicted_field_name, submitted_field_name)

    submitted_genotype_counts = get_genotype_counts(data_dict, field_name=submitted_field_name)
    predicted_genotype_counts = get_genotype_counts(data_dict, field_name=predicted_field_name)

    logging.info("Calcualting F1 for {} samples".format(num_samples))
    raw_scheme_f1 = f1_score(labels['predicted'], labels['submitted'], average='micro')
    num_predicted = num_samples - len(filter_genotypes_exclude(labels,name='')['sample_ids'])
    results = {'overall': {'num_submitted':num_samples,'num_predicted':num_predicted,'f1_score':raw_scheme_f1}}
    for genotype in genotype_rules:
        tmp = filter_genotypes_include(labels, genotype)
        results[genotype] = {    'num_submitted':0,'num_predicted':0,'f1_score':0}
        if len(tmp['sample_ids']) == 0:
            continue
        results[genotype] = {'num_submitted':submitted_genotype_counts[genotype],
                             'num_predicted':predicted_genotype_counts[genotype],
                             'f1_score': f1_score(tmp['predicted'], tmp['submitted'], average='micro')}

    scheme_f1_summary_file = os.path.join(outdir,"{}-scheme.scores.txt".format(prefix))

    write_scheme_scores(SCHEME_SCORES_REPORT_HEADER, scheme_f1_summary_file, os.path.basename(scheme_file), results)

    scheme_sample_file = os.path.join(outdir, "{}-sample.results.txt".format(prefix))
    write_updated_genotype_report(SCHEME_SAMPLE_REPORT_HEADER, scheme_sample_file, os.path.basename(scheme_file),
                                  variants, data_dict, genotype_rules, predicted_field_name,
                                  submitted_field_name)


    logging.info("Analysis complete")

