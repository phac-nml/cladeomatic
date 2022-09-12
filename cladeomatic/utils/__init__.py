import logging
from cladeomatic.constants import LOG_FORMAT
from scipy.stats import entropy, fisher_exact
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
from subprocess import Popen, PIPE
import pandas as pd

def init_console_logger(lvl):
    """

    Parameters
    ----------
    lvl [int] : Integer of level of logging desired 0,1,2,3

    Returns
    -------

    logging object

    """
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[lvl]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging


def run_command(command):
    p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')
    return stdout, stderr

def calc_shanon_entropy(value_list):
    total = sum(value_list)
    if total == 0:
        return -1
    values = []
    for v in value_list:
        values.append(v / total)
    return entropy(values)

def calc_AMI(category_1, category_2):
    return adjusted_mutual_info_score(category_1, category_2, average_method='arithmetic')

def calc_ARI(category_1, category_2):
    return adjusted_rand_score(category_1, category_2)

def parse_metadata(file):
    '''
    Parses metadata file into a dict with samples as the keys
    :param file: str path to metdata file with 'sample_id' as a colum
    :return: dict
    '''
    df = pd.read_csv(file, sep="\t", header=0)
    if not 'sample_id' in df.columns.tolist():
        return {}
    columns = set(df.columns.tolist()) - set('sample_id')
    metadata = {}
    for index, row in df.iterrows():
        metadata[row['sample_id']] = {}
        for field in columns:
            field = str(field)
            if field == 'sample_id':
                continue
            value = str(row[field])
            metadata[row['sample_id']][field] = value
    return metadata