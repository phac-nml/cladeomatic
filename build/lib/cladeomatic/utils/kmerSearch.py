import re
import os
import sys
from itertools import product
import ray
from Bio import SeqIO
from ahocorasick import Automaton

from cladeomatic.utils.seqdata import create_aln_pos_from_unalign_pos_lookup

bases_dict = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'], }

REGEX_GZIPPED = re.compile(r'^.+\.gz$')

NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')

def expand_degenerate_bases(seq):
    """List all possible kmers for a scheme given a degenerate base

    Args:
         Scheme_kmers from SNV scheme fasta file
    Returns:
         List of all possible kmers given a degenerate base or not
    """
    return list(map("".join, product(*map(bases_dict.get, seq))))

def revcomp(s):
    """Reverse complement nucleotide sequence

    Args:
        s (str): nucleotide sequence

    Returns:
        str: reverse complement of `s` nucleotide sequence
    """
    return s.translate(NT_SUB)[::-1]

def init_automaton_dict(seqs):
    """Initialize Aho-Corasick Automaton with kmers from SNV scheme fasta

    Args:
        scheme_fasta: SNV scheme fasta file path

    Returns:
         Aho-Corasick Automaton with kmers loaded
    """
    A = Automaton()
    for seq_id in seqs:
        sequence = seqs[seq_id]
        kmer_list = expand_degenerate_bases(sequence.replace('-',''))
        for idx,seq in enumerate(kmer_list):
            A.add_word(seq, (seq_id, seq, False))
            A.add_word(revcomp(seq), (seq_id, seq, True))

    A.make_automaton()
    return A


@ray.remote
def processSeq(fasta_file, out_file, seqids, num_kmers, klen, aho):
    """
    This method processes the input list of all possible kmers for a temporary
    output file for further processing.  This file contains the sequence id, kmer
    index, kmer start and end indexes, and if the kmer is a reverse compliment.
    Note this method is added to the Ray instance for distributed computational power.
    :param fasta_file: String - the file path to the fasta formatted file.
    :param out_file: String - the path to the output kmer file - temporary processing file
    :param seqids: dictionary - the sequence ids as keys
    :param num_kmers: int - number of kmers to process
    :param klen: int - length of the kmers
    :param aho: Automaton object - the aho-corsick object of the filtered kmer
    """
    #set for reveres complement kmers
    revcomp_kmers = set()
    kmer_align_start = [-1] * num_kmers
    kmer_align_end = [-1] * num_kmers
    seqs_present = {}
    fh = open(out_file,'w')
    #read and process the fasta sequences
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = str(record.id)
            if id not in seqids:
                continue
            #if the sequence id from the fasta is in the
            seq = str(record.seq)
            aln_lookup = create_aln_pos_from_unalign_pos_lookup(seq)
            seq = seq.replace('-', '')
            counts = [0] * num_kmers
            #loop through the Automaton object iterator - from Aho-Corsick
            for idx, (kIndex, kmer_seq, is_revcomp) in aho.iter(seq):
                kIndex = int(kIndex)
                counts[kIndex] += 1
                #if this kmer is a reverse complement, add it to the list
                if is_revcomp:
                    revcomp_kmers.add(kIndex)
                uStart = idx - klen + 1
                uEnd = idx
                #get the start and end indexes for the kmers in this sequence
                if kmer_align_start[kIndex] == -1:
                    kmer_align_start[kIndex] = aln_lookup[uStart]
                    kmer_align_end[kIndex] = aln_lookup[uEnd]
            #loop through the number of kmers to write the temporary file
            #with the sequence id, the kIndex, the count, the kmer start and end positions
            #and if the kmer selected is the reverse compliment
            for kIndex, count in enumerate(counts):
                is_revcomp = kIndex in revcomp_kmers
                row = [
                    id, kIndex, count, kmer_align_start[kIndex], kmer_align_end[kIndex],is_revcomp,
                ]
                fh.write("{}\n".format("\t".join([str(x) for x in row])))
            seqs_present[id] = counts
    fh.close()
    handle.close()


def SeqSearchController(seqKmers, fasta_file,out_dir,prefix,n_threads=1):
    """
    This method takes in the list of all possible kmers for the sequences
    provided and writes them to a temporary processing file for downstream
    searches
    :param seqKmers: dictionary - the index as keys and sequences of the kmers
    :param fasta_file: String - the file path to the fasta sequence file
    :param out_dir: String - the file path to the output file
    :param prefix: String - the prefix for the output files
    :param n_threads: int - number of threads to be used in this process, default is 1
    :return: list - a list of Strings for the paths to the temporary processing file
    """
    num_kmers = len(seqKmers)
    if num_kmers == 0:
        return {}
    klen = len(seqKmers[0])
    count_seqs = 0
    seq_ids = []
    #retrieve the sequence ids from the fasta file
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_ids.append(str(record.id))
            count_seqs+=1
    handle.close()
    #determine the batch size
    batch_size = int(count_seqs / n_threads)
    num_workers = n_threads
    batches = []
    #determine the number of batches based on how many threads and the
    #batch size calculated above
    for i in range(0,num_workers):
        batches.append(seq_ids[i*batch_size:i*batch_size+batch_size])
    #use Ray to find the specific kmers for the from the dictionary of kmers
    aho = ray.put(init_automaton_dict(seqKmers))
    result_ids = []
    file_paths = []
    #loop through each batch to create the temporary kmer file
    for i in range(0,len(batches)):
        outfile = os.path.join(out_dir,"{}-kmersearch-{}.txt".format(prefix,i))
        result_ids.append(processSeq.remote(fasta_file, outfile, batches[i], num_kmers, klen, aho))
        file_paths.append(outfile)
    #Use Ray to set the result ids
    ray.get(result_ids)
    return(file_paths)





