from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
import pandas as pd
from cladeomatic.version import __version__
from cladeomatic.utils import parse_metadata

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
        description="Clade-O-Matic: Genotyping scheme genotype namer v. {}".format(__version__),
        formatter_class=CustomFormatter)
    parser.add_argument('--in_scheme', type=str, required=True,
                        help='Cladeomatic scheme file')
    parser.add_argument('--in_names', type=str, required=True, help='Tab delimited file of (node, name)', default=None)
    parser.add_argument('--outfile', type=str, required=True, help='Output file for updated scheme')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()



def rename(lookup,queries):
    """
    This method takes in the naming data dictionary and the positive or partial genotypes
    list to create a new list of names for the genotypes if they names are the same
    in both the naming dictionary and positive or partial genotypes list.

    Parameters
    ----------
    lookup : dict
        The naming data dictionary
    queries : list
        The list of positive genotypes

    Returns
    -------
    list
        The list of names that are the same between the naming dictionary and the genotypes

    """
    out = []
    for name in queries:
        if name in lookup:
            name = lookup[name]['name']
        out.append(name)
    return out

def run():
    """
    This method takes in the data dictionaries for the scheme and names
    and renames the genotypes accordingly.
    """
    #read the arguments
    cmd_args = parse_args()
    in_scheme = cmd_args.in_scheme
    in_names = cmd_args.in_names
    outfile = cmd_args.outfile
    #read in the naming file
    name_data = parse_metadata(in_names,colum='node')

    #read in the scheme into a dataframe
    df = pd.read_csv(in_scheme, sep="\t", header=0)
    #loop through the rows of the dataframe
    for index, row in df.iterrows():
        #retrieve the positive genotypes from the scheme and split them if they are a list
        positive_genotypes = row['positive_genotypes']
        if isinstance(positive_genotypes,float):
            positive_genotypes = []
        else:
            positive_genotypes = positive_genotypes.split(',')
        #rename the positive genotypes
        positive_genotypes = rename(name_data,positive_genotypes)

        partial_genotypes = row['partial_genotypes']
        if isinstance(partial_genotypes,float):
            partial_genotypes = []
        else:
            partial_genotypes = partial_genotypes.split(',')
        partial_genotypes = rename(name_data, partial_genotypes)
        row['positive_genotypes'] = ",".join([str(x) for x in positive_genotypes])
        row['partial_genotypes'] = ",".join([str(x) for x in partial_genotypes])
        df.loc[index, ['positive_genotypes','partial_genotypes']] = [row['positive_genotypes'],row['partial_genotypes']]

    df.to_csv(outfile,sep="\t",header=True,index=False)

