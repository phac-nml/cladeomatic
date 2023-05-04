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
from deprecated import deprecated


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
                    Instantiates the default values for the ArgumentParser for display on the command line.
                RawDescriptionHelpFormatter object
                    Ensures the correct display of the default values for the ArgumentParser
                """
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
    """
    This method retrieves both the SNP and sample profiles from the
    VCF file, ensures the SNPs are in valid positions and adds the sample
    SNP profile to the data dictionary.

    Parameters
    ----------
    valid_positions : list
        The list of integers indicating the valid SNP positions for this sample set
    vcf_file : str
        The file path to vcf or tsv snp data files

    Returns
    -------
    dict
        The dictionary of the snp data
    """
    vcf = vcfReader(vcf_file)
    #read the first row of the VCF file
    data = vcf.process_row()
    #retrieve the list of samples from the vcf file
    samples = vcf.samples
    profiles = {}
    count_snps = 0
    if data is not None:
        count_snps += 1
    #initialize the sample snp profiles data dictionary
    for sample_id in samples:
        profiles[sample_id] = {}
    if data is None:
        return profiles
    #loop through the vcf data file rows
    while data is not None:
        pos = int(data['POS'])
        #ensure the position of the snp from the VCF file is in the valid
        #positions list
        if not pos in valid_positions:
            data = vcf.process_row()
            continue
        #if the snp is valid, add the snp to the sample snp profile data dictionary
        for sample_id in samples:
            base = data[sample_id]
            profiles[sample_id][pos] = base
        count_snps += 1

        data = vcf.process_row()

    return profiles

def parse_scheme_features(scheme_file):
    """
    This method parses the pass snp scheme file and creates a scheme and feature data set for
    downstream processing.  The scheme data dictionary consists of the snp scheme data
    of the snp position, the snp base and if the genotype is positive, partial and/or allowed.
    The feature data dictionary the number of positions, the number of mutations, total features
    and the mutation data for searching.  This feature data dictionary also included the features
    of each genotype: the mutation lookup id, if the snp was a ref or alt state, the base
    and base position.

    Parameters
    ----------
    scheme_file : str
        The file path to the snp scheme file to read in

    Returns
    -------
    dict
        The constructed features data dictionary, with the scheme data dictionary inside

    """
    features = {
        'scheme':{},
        'num_positions':0,
        'num_mutations':0,
        'total_scheme_features':0,
        'mutation_lookup':{}
    }

    scheme = {}
    #read the scheme file into a pandas dataframe
    df = pd.read_csv(scheme_file, sep="\t", header=0, low_memory=False)
    snp_states = {}
    #get the various features and set them in the features dictionary
    num_positions = df['variant_start'].nunique()
    num_features = df['mutation_key'].nunique()
    features['total_scheme_features'] = len(df)
    features['num_positions'] = num_positions
    features['num_mutations'] = num_features

    #loop through the rows of the dataframe
    for row in df.itertuples():
        #get the state (ref or alt)
        state = row.state
        mutation_key = row.mutation_key
        target_variant = row.target_variant
        variant_start = int(row.variant_start)
        if not variant_start in features['mutation_lookup']:
            features['mutation_lookup'][variant_start] = {'alt':{},'ref':{}}
        features['mutation_lookup'][variant_start][state][target_variant] = mutation_key

        positive_genotypes = row.positive_genotypes
        #check to see if the positive genotypes are a single float for one sample or a list of sample genotypes
        #if a list, then split the list into a list
        if isinstance(positive_genotypes, float):
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')

        genotypes = positive_genotypes
        #loop through the list of positive genotypes
        for genotype in positive_genotypes:
            #add the genotype to the scheme if it hasn't already been added
            if not genotype in scheme:
                scheme[genotype] = {'positive': {}, 'partial': {}, 'allowed': {}}
            #if the start of the variant for this genotype has not been added, add it
            if not variant_start in scheme[genotype]['positive']:
                scheme[genotype]['positive'][variant_start] = set()
            scheme[genotype]['positive'][variant_start].add(target_variant)
        #get the partial genotypes
        partial_genotypes = row.partial_genotypes
        #add the variant start position to the snp state
        if not variant_start in snp_states:
            snp_states[variant_start] = set()
        snp_states[variant_start].add(target_variant)
        # check to see if the genotypes are a single float for one sample or a list of sample genotypes
        # if a list, then split the list into a list
        if isinstance(partial_genotypes, float):
            partial_genotypes = []
        else:
            partial_genotypes = partial_genotypes.split(',')
        genotypes += partial_genotypes
        if len(partial_genotypes) == 0:
            continue
        #loop through partial genotypes to add them to the scheme data set
        for genotype in partial_genotypes:
            if not genotype in scheme:
                scheme[genotype] = {'positive': {}, 'partial': {}, 'allowed': {}}
            if not variant_start in scheme[genotype]['partial']:
                scheme[genotype]['partial'][variant_start] = set()
            scheme[genotype]['partial'][variant_start].add(target_variant)
    #loop through genotypes added to the scheme data set to
    #add the allowed bases for the genotype position
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


@deprecated
def parse_scheme_genotypes(scheme_file):
    """
    This method parses the snp scheme file to construct a scheme data dictionary
    for further downstream processing.

    Parameters
    ----------
    scheme_file : str
        The file path the snp scheme file for parsing

    Returns
    -------
    dict
        The scheme data dictionary

    """
    scheme = {}
    df = pd.read_csv(scheme_file, sep="\t", header=0, low_memory=False)
    snp_states = {}
    #read the snp scheme data frame
    for row in df.itertuples():
        target_variant = row.target_variant
        variant_start = int(row.variant_start)
        positive_genotypes = row.positive_genotypes
        if isinstance(positive_genotypes,float):
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')
        #add the genotypes to the scheme
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
        #add the genotypes to the scheme
        for genotype in partial_genotypes:
            if not genotype in scheme:
                scheme[genotype] = {'positive':{},'partial':{},'allowed':{}}
            if not variant_start in scheme[genotype]['partial']:
                scheme[genotype]['partial'][variant_start] = set()
            scheme[genotype]['partial'][variant_start].add(target_variant)

    #loop through the genotypes in the scheme to add the allowed based for the genotype variants
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
    """
    This method is used to call the genotypes of the created scheme.
    These calls are based on the snp scheme data created in the
    :meth:`parse_scheme_features` method.  This method determines
    which genotypes are valid based on the scheme data dictionary
    and the variants

    Parameters
    ----------
    genotype_rules : dict
        The snp scheme data dictionary
    metadata : dict
        The sample metadata dictionary
    variants : dict
        The snp variants data dictionary
    max_dist : int
        The maximum distance between alleles

    Returns
    -------
    dict
        The results of the genotyping call in the form of a dictionary.  This dictionary consists of the genotype, its allelic distance, and the matching snps.
    Notes
    -----
    Refer to https://www.ray.io for more information about the Ray instances used in this module.
    """
    result = {}
    #loop through the sample metadata
    for sample_id in metadata:
        #if the smple is not the variants list, skip
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
        #initialize the genotype results dictionary with an entry for each
        #genotype entry in the rules dictionary
        for genotype in genotype_rules:
            genoytpe_results[genotype] = {'match':{},'mismatch':{}}
            dists[genotype] = 1
        #loop through the snp positions in the variant sample dictionary
        for pos in variants[sample_id]:
            found_base = set(variants[sample_id][pos])
            #skip invalid bases
            if '*' in found_base  or '-' in found_base or 'N' in found_base:
                continue
            #loop through the genotypes
            for genotype in genotype_rules:
                #if the allowed bases in the genotype
                allowed_bases = genotype_rules[genotype]['allowed'][pos]
                #if the bases are the same, this is a match, else it is a mismatch
                #set the results in the genotype results data dictionary
                if len(found_base & allowed_bases) == 1:
                    genoytpe_results[genotype]['match'][pos] = found_base

                else:
                    genoytpe_results[genotype]['mismatch'][pos] = found_base

        #loop through the genotype results
        #to calculate the distance between the matching and mismatching
        #alleles
        for genotype in genoytpe_results:
            matches = len(genoytpe_results[genotype]['match'])
            mismatches = len(genoytpe_results[genotype]['mismatch'])

            total = matches + mismatches
            if total > 0:
                dists[genotype] = 1 - matches /total

        #sorts the distances dictionary and keeps the same keys (sample ids)
        result[sample_id]['genoytpe_dists'] =  {k: v for k, v in sorted(dists.items(), key=lambda item: item[1])}

        pdist = 1

        #loop through the genotypes distances
        for genotype in result[sample_id]['genoytpe_dists']:
            dist =  result[sample_id]['genoytpe_dists'][genotype]
            #append only those results that are smaller than the pdist
            #which is updated every allele iteration
            #and the max_dist
            if dist <= pdist and dist <= max_dist:
                result[sample_id]['predicted_genotype(s)'].append(genotype)
                result[sample_id]['predicted_genotype_dist'].append(dist)
                result[sample_id]['genoytpe_results'].append({'match':genoytpe_results[genotype]['match'],
                                    'mismatch':genoytpe_results[genotype]['mismatch']})
                pdist = dist
        #delete the genotype distances as we have the distances we are looking for
        del result[sample_id]['genoytpe_dists']
    return result


def convert_features_to_mutations(feature_lookup,snp_profile):
    """
    This method converts the mutations or features in the scheme
    to the snp base mutation and sets this on a data dictionary
    with the snp state.

    Parameters
    ----------
    feature_lookup : dict
        The mutations or features in the scheme data dictionary
    snp_profile : dict
        The snp profile for a given genotype

    Returns
    -------
    dict
        A dictionary of the snp ids and if they belong to the alt or ref for the genotype
    """
    found_features = {
        'alt':set(),
        'ref':set()
    }
    #loop through the snp profiles for the snp positions
    for pos in snp_profile:
        #skip any entries not in the mutation map
        if pos not in feature_lookup:
            continue
        base = snp_profile[pos]
        #loop through the states for the position
        for state in feature_lookup[pos]:
            #if the base from the snp profile is found in the mutation profile for the state
            #add the base to the found features dictionary for the state
            if base in feature_lookup[pos][state]:
                found_features[state].add(feature_lookup[pos][state][base])
                break

    return found_features


def write_genotype_calls(header,scheme_name,outfile,genotype_results,sample_metadata,genotype_meta, scheme_data, sample_variants,min_positions=1):
    """
    This method writes the genotype calls determined through previous methods to
    a file.  Please refer to the file
    examples/small_test/cladeomatic/cladeomatic-genotype.calls.txt for more information.

    Parameters
    ----------
    header : str
        The header for the output file
    scheme_name : str
        The name for the scheme
    outfile : str
        The output file path
    genotype_results : dict
        The dictionary of the genotype calls from :meth:`call_genotypes`
    sample_metadata : dict
        The sample metadata data dictionary
    genotype_meta : dict
        The genotype metadata data dictionary
    scheme_data : dict
        The snp scheme data dictionary
    sample_variants : dict
        The sample variants data dictionary
    min_positions : int
        The minimum number of . Default of 1.

    """
    #get all additional_fields for writing
    sample_fields = set()
    #get the sample metadata fields
    for sample_id in sample_metadata:
        for field in sample_metadata[sample_id]:
            sample_fields.add(field)
    sample_fields = sorted(list(sample_fields))
    genotype_fields = set()
    #get the genotype metadata fields
    for genotype in genotype_meta:
        for field in genotype_meta[genotype]:
            genotype_fields.add(field)
    genotype_fields = sorted(list(genotype_fields))
    #combine the header, sample and genotype metadata fields
    header = header + sample_fields + genotype_fields

    #initialize file
    fh = open(outfile, 'w')
    fh.write("{}\n".format("\t".join([str(x) for x in header])))
    #create the analysis date
    analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    num_positions = scheme_data['num_positions']
    total_scheme_features = scheme_data['total_scheme_features']

    #loop through the genotype call results and add the data to a row dictionary
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
        #determine the total number of positions by finding the total number of unique mutations in both the ref and alt states
        row['num_detected_positions'] =  len(detected_mutations['alt'] | detected_mutations['ref'])
        #set the number of alt mutations
        row['detected_alt_mutations'] = ";".join([str(x) for x in sorted(list(detected_mutations['alt']))])

        #set the filed metadata to the row
        for field_id in sample_metadata[sample_id]:
            row[field_id] = sample_metadata[sample_id][field_id]

        #if there is only one predicted genotype for the sample id, then set the
        #predicted genotype, its distance and metadata to the row and the genotype
        if len(genotype_results[sample_id]['predicted_genotype(s)']) == 1:
            row['predicted_genotype'] = genotype_results[sample_id]['predicted_genotype(s)'][0]
            row['predicted_genotype_distance'] = genotype_results[sample_id]['predicted_genotype_dist'][0]
            genotype = str(row['predicted_genotype'])
            if genotype in genotype_meta:
                for field_id in genotype_meta[genotype]:
                    row[field_id] = genotype_meta[genotype][field_id]
        #if there are more than one predicted genotypes, set a quality control message
        elif len(genotype_results[sample_id]['predicted_genotype(s)']) > 1:
            status = 'Warning'
            row['qc_messages'] = "Ambiguous genotype assignment, possible genotypes: {}".format(";".join([str(x) for x in genotype_results[sample_id]['predicted_genotype(s)']]))
        else:
            status = 'Warning'
            row['qc_messages'] = 'No genotypes were compatible with the sample'

        #if there are fewer detected positions for snps than the threshold for validity,
        #set an error message in the call file
        if row['num_detected_positions'] < min_positions:
            status = 'Fail'
            row['qc_messages'] = 'Sample is missing too many positions for genotype call'

        row['qc_status'] = status
        fh.write("{}\n".format("\t".join([str(x) for x in row.values()])))

    fh.close()

def is_valid_file(filename):
    """
    A helper method to determine if the file path leads to a valid file.  For the
    file to be considered valid it must exist and have the minimum file size.

    Parameters
    ----------
    filename:str
        The path to the file requiring validation

    Returns
    -------
    bool
        True if this is a valid file, False if not
    """
    status = True
    if not os.path.isfile(filename) or os.path.getsize(filename) <= MIN_FILE_SIZE :
        status = False
    return status


def run():
    """
    The main method to read the command line arguments and creates the genotype call file.
    This method reads the scheme, variant, sample and genotype metadata files, and constructs
    various processing data dictionaries to ultimately call the genotypes and provide a measure of
    quality control of those calls for the samples passed.

    Notes
    -----
    Refer to https://www.ray.io for more information about the Ray instances used in this module.

    """
    #read in the arguments
    cmd_args = parse_args()
    scheme_file = cmd_args.in_scheme
    variant_file = cmd_args.in_var
    metadata_file = cmd_args.sample_meta
    genotype_meta_file = cmd_args.genotype_meta
    outfile = cmd_args.outfile
    num_threads = cmd_args.num_threads
    max_missing_positions = cmd_args.max_missing_positions
    #initialize ray
    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True, num_cpus=num_threads)

    #initialize the logger
    logging = init_console_logger(3)
    logging.info("Starting analysis")

    #validate and read the sample metadata
    logging.info("Reading metadata file {}".format(metadata_file))
    if not is_valid_file(metadata_file):
        logging.error("Error file {} was not found or is empty".format(metadata_file))
        sys.exit()
    sample_metadata = parse_metadata(metadata_file)

    #validate and read the genotype metadata
    genotype_metadata = {}
    if genotype_meta_file is not None:
        logging.info("Reading metadata file {}".format(genotype_meta_file))
        if not is_valid_file(genotype_meta_file):
            logging.error("Error file {} was not found or is empty".format(genotype_meta_file))
            sys.exit()
        logging.info("Reading metadata file {}".format(genotype_meta_file))
        genotype_metadata = parse_metadata(genotype_meta_file,column='key')

    #validate and read the scheme data
    logging.info("Reading scheme file {}".format(scheme_file))
    if not is_valid_file(scheme_file):
        logging.error("Error file {} was not found or is empty".format(scheme_file))
        sys.exit()

    scheme_data = parse_scheme_features(scheme_file)
    genotype_rules = scheme_data['scheme']
    #validate the scheme data
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

    #call the genotypes using ray
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
    #write the genotype calls
    write_genotype_calls(GENOTYPE_REPORT_HEADER, os.path.basename(scheme_file), outfile, genotype_results, sample_metadata, genotype_metadata, scheme_data,variants, min_positions)

    logging.info("Analysis complete")

