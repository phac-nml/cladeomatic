from subprocess import Popen, PIPE

def run_snpdists(aln_file,out_file,num_threads=1):
    """
    Runs the command line snp-dists application which calculates the SNP
    distances from a fasta file.

    Parameters
    ----------
    aln_file : str
        The path to the fasta sequence file
    out_file : str
        The path to the output file
    num_threads : int
        The number of processing threads for . Default is 1.

    Returns
    -------
    tuple
        Stderr - the console error messages

    Notes
    -----
    Please refer to https://github.com/tseemann/snp-dists for more documentation
    """
    command = [
        'snp-dists',
        '-j',num_threads,
        aln_file,
    ]
    fh = open(out_file,'wb')
    cmd = " ".join([str(x) for x in command])
    p = Popen(cmd, shell=True, stdout=fh, stderr=PIPE)
    stderr = p.communicate()
    fh.close()
    return stderr
