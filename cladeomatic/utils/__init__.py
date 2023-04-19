import logging
from cladeomatic.constants import LOG_FORMAT
from scipy.stats import entropy, fisher_exact
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
from subprocess import Popen, PIPE
import pandas as pd
from deprecated import deprecated

def init_console_logger(lvl):
    """
    This method initializes the logger for the console
    :param lvl: int - level of logging desired 0,1,2,3
    :return: logging - the logging object
    """
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[lvl]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging


def run_command(command):
    """
    This method runs the passed command on the shell command line
    :param command: String - the command for the command line
    :return: bytes - stdout, stderr: the standard out and error messages
    returned by the command line
    """
    p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')
    return stdout, stderr

def calc_shanon_entropy(value_list):
    """
    This method calculates the shannon entropy value for the list
    of numbers passed
    :param value_list: list - the list of values to use for the calculation
    :return: float - the calculated Shannon entropy or -1 if there are no
    values in the passed list
    """
    total = sum(value_list)
    if total == 0:
        return -1
    values = []
    for v in value_list:
        values.append(v / total)
    return entropy(values)
@deprecated()
def calc_AMI(category_1, category_2):
    """
    Calculates the adjusted mutual info score between two clusterings
    :param category_1: list - a list of int values for cluster 1
    :param category_2: list - a list of int values for cluster 2
    :return: float - the value of the AMI score for the two clusters
    """
    return adjusted_mutual_info_score(category_1, category_2, average_method='arithmetic')
@deprecated
def calc_ARI(category_1, category_2):
    """
    Calculates the adjusted rand score between two clusterings
    :param category_1: list - a list of int values for cluster 1
    :param category_2: list - a list of int values for cluster 2
    :return: float - the value of the ARI score for the two clusters
    """
    return adjusted_rand_score(category_1, category_2)

def parse_metadata(file):
    """
    Parses the metadata file into a dictionary with sample ids as the keys.
    Will dynamically add the other columns and values to the dictionary
    :param file: string - path to tsv metadata file with 'sample_id' as a
    mandatory column
    :return: dictionary - a dictionary of the parsed metadata values organized
    wtih the sample id as a key
    """
    #read the file into a pandas datafram
    df = pd.read_csv(file, sep="\t", header=0)
    if not 'sample_id' in df.columns.tolist():
        return {}
    columns = set(df.columns.tolist()) - set('sample_id')
    metadata = {}
    for index, row in df.iterrows():
        metadata[row['sample_id']] = {}
        #parse out any additional columns and their row values
        for field in columns:
            field = str(field)
            if field == 'sample_id':
                continue
            value = str(row[field])
            metadata[row['sample_id']][field] = value
    return metadata