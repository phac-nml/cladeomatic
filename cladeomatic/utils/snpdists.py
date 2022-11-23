from subprocess import Popen, PIPE

def run_snpdists(aln_file,out_file,num_threads=1):
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
