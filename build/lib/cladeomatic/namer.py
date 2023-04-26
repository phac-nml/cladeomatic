from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
import pandas as pd
from cladeomatic.version import __version__


def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Clade-O-Matic: Genotyping scheme genotype namer v. {}".format(__version__),
        formatter_class=CustomFormatter)
    parser.add_argument('--in_scheme', type=str, required=True,
                        help='Cladeomatic scheme file')
    parser.add_argument('--in_names', type=str, required=True, help='Tab delimited file of (genotype, name)', default=None)
    parser.add_argument('--outfile', type=str, required=True, help='Output file for updated scheme')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()

def parse_names(file):
    '''
    Parses metadata file into a dict with genotype as the keys
    :param file: str path to metdata file with 'genotype and name' as a columns
    :return: dict
    '''
    df = pd.read_csv(file, sep="\t", header=0)
    if not 'genotype' in df.columns.tolist():
        return {}
    metadata = {}
    for index, row in df.iterrows():
        metadata[row['genotype']] = row['name']
    return metadata

def rename(lookup,queries):
    out = []
    for name in queries:
        if name in lookup:
            name = lookup[name]
        out.append(name)
    return out

def main():
    cmd_args = parse_args()
    in_scheme = cmd_args.in_scheme
    in_names = cmd_args.in_names
    outfile = cmd_args.outfile

    names = parse_names(in_names)
    df = pd.read_csv(in_scheme, sep="\t", header=0)
    for index, row in df.iterrows():
        positive_genotypes = row['positive_genotypes']
        if isinstance(positive_genotypes,float):
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')
        positive_genotypes = rename(names,positive_genotypes)

        partial_genotypes = row['partial_genotypes']
        if isinstance(partial_genotypes,float):
            partial_genotypes = []
        else:
            partial_genotypes = partial_genotypes.split(',')
        partial_genotypes = rename(names, partial_genotypes)
        row['positive_genotypes'] = ",".join([str(x) for x in positive_genotypes])
        row['partial_genotypes'] = ",".join([str(x) for x in partial_genotypes])
        df.loc[index, ['positive_genotypes','partial_genotypes']] = [row['positive_genotypes'],row['partial_genotypes']]

    df.to_csv(outfile,sep="\t",header=True,index=False)


    return