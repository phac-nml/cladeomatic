from subprocess import Popen, PIPE

def run_snpdists(aln_file,out_file,num_threads=1):
    """
    Runs the command line snp-dists application which calculates the SNP
    distances from a fasta file.
    :param aln_file: String - the path to the fasta sequence file
    :param out_file: String - the path to the output file
    :param num_threads: int - number of processing threads, default is 1
    :return: tuple - Stderr - the console error messages
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
