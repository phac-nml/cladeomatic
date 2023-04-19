import pandas as pd
from cladeomatic.utils import run_command
import pandas as pd
from deprecated import deprecated

from cladeomatic.utils import run_command

@deprecated
def run_jellyfish_count(seq_file,out_file,mem='10M',k=21,n_threads=1,jf_out=False):
    if jf_out:
        cmd = "jellyfish count -t {} -m {} -s {} -C -o {} {}".format(n_threads, k, mem, out_file, seq_file)
    else:
        cmd = "jellyfish count --text -t {} -m {} -s {} -C -o {} {}".format(n_threads,k,mem,out_file,seq_file)

    (stdout,stderr) = run_command(cmd)
    print(stdout)
    print(stderr)
    return (stdout,stderr)
@deprecated()
def run_jellyfish_query(seq_file,jellyfish_mers,out_file):
    cmd = "jellyfish query -o {} {} {}".format(out_file, seq_file, jellyfish_mers)
    (stdout,stderr) = run_command(cmd)
    return (stdout,stderr)


@deprecated()
def parse_jellyfish_counts(file):
    return pd.read_csv(file,skiprows=1, header=None,names=['kmer','count'],sep=' ')


