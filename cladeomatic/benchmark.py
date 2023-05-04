import logging
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
    """ Argument Parsing method.

        A function to parse the command line arguments passed at initialization of Clade-o-matic,
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
    """
    A helper method to read in the genotype call file for downstream use.

    Parameters
    ----------
    header : list
        The list of strings for the genotype call report header
    file : str
        The path to the genotype call report file to parse

    Returns
    -------
    dict
        A dictionary of the parsed genotype call report file
    """
    return parse_metadata(file,column=header[0])

def convert_genotype_report(data_dict,predicted_field_name,submitted_field_name):
    """
    A method to take in the genotype call data dictionary, the predicted
    and submitted field names and converts these to a data dictionary
    consisting of the sample ids for the genotypes, the predicted genotype,
    and submitted genotype ids.

    Parameters
    ----------
    data_dict : dictionary
        The genotype call dictionary to parse
    predicted_field_name : str
        The name of the predicted genotype field within the genotype call file
    submitted_field_name :str
        The name of the submitted genotype field within the genotype call file

    Returns
    -------
    dict
        A dictionary of sample_ids, the predicted genotype calls and the submitted genotype calls
    """
    p = []
    s = []
    samples = []
    #loop through the genotype call data dictionary and create the data dictionary
    for sample_id in data_dict:
        samples.append(str(sample_id))
        p.append(str(data_dict[sample_id][predicted_field_name]))
        s.append(str(data_dict[sample_id][submitted_field_name]))
    return {'sample_ids':samples,'predicted':p,'submitted':s}

def filter_genotypes_include(labels,name):
    """
    This method takes the dictionary of sample_ids,
    the predicted genotype calls and the submitted genotype calls,
    and filters the dictionary for the genotype name.  If the genotype
    name exists in either the predicted or the submitted genotype call
    the sample will be added to the dictionary and the inclusive
    dictionary is returned.

    Parameters
    ----------
    labels : dict
        A dictionary of sample_ids, the predicted genotype calls and the submitted genotype calls
    name : str
        The name of the genotype to filter by

    Returns
    -------
    dict
        The genotype name-filtered dictionary of sample_ids, the predicted genotype calls and the submitted genotype calls

    """
    p = []
    s = []
    samples = []
    #loop through the data dictionary
    for i in range(0,len(labels['sample_ids'])):
        #if the genotype name is in either the submitted or predicted gentoype call
        #columns, then add the sample id, the predicted and submitted genotypes
        #to the data dictionary
        if name in labels['predicted'][i] or name in labels['submitted'][i]:
            samples.append(labels['sample_ids'][i])
            p.append(labels['predicted'][i])
            s.append(labels['submitted'][i])
    return {'sample_ids':samples,'predicted':p,'submitted':s}

def filter_genotypes_exclude(labels,name):
    """
    This method takes the dictionary of sample_ids,
    the predicted genotype calls and the submitted genotype calls,
    and filters the dictionary for the genotype name.  If the genotype
    name exists in either the predicted or the submitted genotype call
    the sample is excluded from addition to the data dictionary.

    Parameters
    ----------
    labels : dict
        A dictionary of sample_ids, the predicted genotype calls and the submitted genotype calls
    name : str
        The name of the genotype to filter by

    Returns
    -------
    dict
        The genotype name-filtered dictionary of sample_ids, the predicted genotype calls and the submitted genotype calls

    """
    p = []
    s = []
    samples = []
    # loop through the data dictionary
    for i in range(0,len(labels['sample_ids'])):
        # if the genotype name is in either the submitted or predicted genetype call
        # columns, skip.
        if name in labels['predicted'][i] or name in labels['submitted'][i]:
            continue
        #only add the samples that do not contain the genotype name
        samples.append(labels['sample_ids'][i])
        p.append(labels['predicted'][i])
        s.append(labels['submitted'][i])
    return {'sample_ids':samples,'predicted':p,'submitted':s}

def write_scheme_scores(header,file,scheme_name,scores):
    """
    This method writes the scheme scores determined to a file.  Please
    refer to the file examples/small_test/cladeomatic/benchmark/cladeomatic-scheme.scores.txt for
    more information.

    Parameters
    ----------
    header : str
        The header string for the output file
    file : str
        The path to the output file
    scheme_name : str
        The name of the scheme
    scores : dict
        The data dictionary of the scores to be written out

    """
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
    """
    This method identifies which snps exist in the variants data dictionary
    but not in the genotype rules or schema dictionary for a specific sample
    id.

    Parameters
    ----------
    profile : dict
        The variants data dictionary for a specific sample id
    rule : dict
        The scheme data dictionary for a specific sample id

    Returns
    -------
    dict
        A data dictionary of the snp positions and bases that are not found in both the variants and scheme dictionaries for the sample passed
    """
    mismatches = {}
    for pos in profile:
        #retireve the snp base
        found_base = set(profile[pos])
        if '*' in found_base or '-' in found_base or 'N' in found_base:
            continue
        #find the allowed bases in the rules dictionary
        allowed_bases = rule['allowed'][pos]
        #if there is a mis-match, add the position and mismatched base
        if len(found_base & allowed_bases) != 1:
            mismatches[pos] = found_base
    return mismatches

def write_updated_genotype_report(header,file,scheme_name,variants,data_dict,genotype_rules,predicted_field_name,submitted_field_name):
    """
    This method writes the genotype sample results file with the sample_id, scheme name,
    submitted genotype, predicted genotype, if the genotypes are a match for each other,
    and if not what are the problem snp positions and bases (features). Please refer to
    the file examples/small_test/cladeomatic/benchmark/cladeomatic-sample.results.txt for
    more information.

    Parameters
    ----------
    header : str
        The header string for the output file
    file : str
        The file path for the output file
    scheme_name : str
        The name for the scheme
    variants : dict
        The data dictionary for the variant call file
    data_dict : dict
        The data dictionary for the called genotypes
    genotype_rules : dict
        The data dictionary for the genotype scheme/rules
    predicted_field_name : str
        The name for the predicted genotype field name
    submitted_field_name : str
        The name for the submitted genotype field name

    """
    fh = open(file,'w')
    fh.write("{}\n".format("\t".join(header)))
    analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    for sample_id in data_dict:
        p = data_dict[sample_id][predicted_field_name]
        s = data_dict[sample_id][submitted_field_name]
        is_match = p == s
        mismatches = {}
        #find the mismatches, which genotype has the problem

        #there is an issue in this code similar to that in the run() method
        #where there is a key error, quick hack around
        if not is_match:
            try:
                mismatches = get_problem_bases(variants[sample_id],genotype_rules[s])
            except KeyError as ke:
                logging.error("Key error {}".format(ke))
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
    """
    This method counts the number of unique genotypes that appear in the data
    dictionary passed.

    Parameters
    ----------
    data_dict : dict
        The data dictionary that contains genotypes for counting
    field_name : str
        The name of the field within the data dictionary for which to parse and count genotypes for.  Default is 'genotype'.

    Returns
    -------
    dict
        The dictionary of genotypes and their counts
    """
    counts = {}
    for sample_id in data_dict:
        genotype = data_dict[sample_id][field_name]
        if not genotype in counts:
            counts[genotype] = 0
        counts[genotype]+=1
    return counts



def run():
    """
    The main method to read the command line arguments and creates a file for the
    F1 scores and if required a file that records per sample any sites
    that are responsible for the genotype not being called.
    This method reads the genotype call file, the variant call file and the
    kmer inclusive scheme file to verify the genotypes and determine if the
    submitted and predicted genotypes are a match.  The results are written
    to the output files
    examples/small_test/cladeomatic/benchmark/cladeomatic-scheme.scores.txt and
    examples/small_test/cladeomatic/benchmark/cladeomatic-sample.results.txt.

    """
    #read in the arguments
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
    #read in the genotype call file
    scheme_data = parse_scheme_features(scheme_file)
    genotype_rules = scheme_data['scheme']
    #loop through the scheme data, if there are no positive kmers, issue a warning via the log
    for genotype in genotype_rules:
        if len(genotype_rules[genotype]['positive']) == 0:
            logging.warn("Genotype {} has no required kmers".format(genotype))

    valid_positions = []
    #loop through the genotypes and aggregate the valid positions (positivem partial and allowed)
    for genotype in genotype_rules:
        for state in genotype_rules[genotype]:
            valid_positions += list(genotype_rules[genotype][state].keys())
    valid_positions = list(set(valid_positions))

    logging.info("Extracted {} genotyping positions".format(len(valid_positions)))

    logging.info("Reading snp data from vcf file {}".format(variant_file))
    variants = get_snp_profiles(valid_positions, variant_file)

    #read in the genotype calls file and simplified it for id, predicted and submitted fields
    data_dict = parse_genotype_report(GENOTYPE_REPORT_HEADER,geno_file)
    num_samples = len(data_dict)
    labels = convert_genotype_report(data_dict, predicted_field_name, submitted_field_name)

    #retrieve the submitted and predicted gene counts from the above simplified genotype call dictionary
    submitted_genotype_counts = get_genotype_counts(data_dict, field_name=submitted_field_name)
    predicted_genotype_counts = get_genotype_counts(data_dict, field_name=predicted_field_name)

    logging.info("Calculating F1 for {} samples".format(num_samples))
    raw_scheme_f1 = f1_score(labels['predicted'], labels['submitted'], average='micro')
    num_predicted = num_samples - len(filter_genotypes_exclude(labels,name='')['sample_ids'])
    results = {'overall': {'num_submitted':num_samples,'num_predicted':num_predicted,'f1_score':raw_scheme_f1}}
    for genotype in genotype_rules:
        tmp = filter_genotypes_include(labels, genotype)
        results[genotype] = {'num_submitted':0,'num_predicted':0,'f1_score':0}
        if len(tmp['sample_ids']) == 0:
            continue
        #there is a known bug in the below code, adding a try catch for now
        try:
            results[genotype] = {'num_submitted':submitted_genotype_counts[genotype],
                                 'num_predicted':predicted_genotype_counts[genotype],
                                 'f1_score': f1_score(tmp['predicted'], tmp['submitted'], average='micro')}
        except KeyError as ke:
            logging.error("Key error on genotype lookup".format(ke))

    scheme_f1_summary_file = os.path.join(outdir,"{}-scheme.scores.txt".format(prefix))
    # write the f1 scores summary file
    write_scheme_scores(SCHEME_SCORES_REPORT_HEADER, scheme_f1_summary_file, os.path.basename(scheme_file), results)

    scheme_sample_file = os.path.join(outdir, "{}-sample.results.txt".format(prefix))
    #write the genotype report
    write_updated_genotype_report(SCHEME_SAMPLE_REPORT_HEADER, scheme_sample_file, os.path.basename(scheme_file),
                                  variants, data_dict, genotype_rules, predicted_field_name,
                                  submitted_field_name)

    logging.info("Analysis complete")

