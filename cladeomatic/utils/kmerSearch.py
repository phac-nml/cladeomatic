import re
from itertools import product
import tempfile
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
def processSeq(fasta_file, seqids, num_kmers, klen, aho):
    invalid_kmers = set()
    revcomp_kmers = set()
    kmer_align_start = [-1] * num_kmers
    kmer_align_end = [-1] * num_kmers
    seqs_present = {}
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = str(record.id)
            if id not in seqids:
                continue
            seq = str(record.seq)
            aln_lookup = create_aln_pos_from_unalign_pos_lookup(seq)
            seq = seq.replace('-', '')
            counts = [0] * num_kmers
            for idx, (kIndex, kmer_seq, is_revcomp) in aho.iter(seq):
                kIndex = int(kIndex)
                counts[kIndex] += 1
                if is_revcomp:
                    revcomp_kmers.add(kIndex)
                uStart = idx - klen + 1
                uEnd = idx
                kmer_align_start[kIndex] = aln_lookup[uStart]
                kmer_align_end[kIndex] = aln_lookup[uEnd]

            for kIndex, count in enumerate(counts):
                if count > 1:
                    invalid_kmers.add(kIndex)
            seqs_present[id] = counts
    handle.close()
    return {
        'invalid_kmers': invalid_kmers,
        'revcomp_kmers': revcomp_kmers,
        'kmer_align_start': kmer_align_start,
        'kmer_align_end': kmer_align_end,
        'seqs_kcounts':seqs_present
    }


def SeqSearchController(seqKmers, fasta_file,out_kmer_file,n_threads=1):
    num_kmers = len(seqKmers)
    if num_kmers == 0:
        return {}
    klen = len(seqKmers[0])
    count_seqs = 0
    seq_ids = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_ids.append(str(record.id))
            count_seqs+=1
    handle.close()
    batch_size = int(count_seqs / n_threads)
    num_workers = n_threads
    batches = []
    for i in range(0,num_workers):
        batches.append(seq_ids[i*batch_size:i*batch_size+batch_size])
    aho = ray.put(init_automaton_dict(seqKmers))
    result_ids = []
    for i in range(0,len(batches)):
        result_ids.append(processSeq.remote(fasta_file, batches[i], num_kmers, klen, aho))

    results = ray.get(result_ids)
    del(aho)
    invalid_kmers = set()
    revcomp_kmers = set()
    kmer_align_start = results[0]['kmer_align_start']
    kmer_align_end = results[0]['kmer_align_end']
    kmer_counts = {}
    for i in range(0,len(results)):
        invalid_kmers = invalid_kmers | results[i]['invalid_kmers']
        revcomp_kmers = revcomp_kmers | results[i]['revcomp_kmers']
        kmer_counts.update(results[i]['seqs_kcounts'])

    filt = {}
    id = 0
    fh = open(out_kmer_file,'w')
    fh.write("id\taStart\taEnd\tkseq\tseqs_present\n")
    for idx in range(0,num_kmers):
        if idx in invalid_kmers:
            continue
        s = kmer_align_start[idx]
        e = kmer_align_end[idx]
        if s == -1 or e == -1:
            continue
        kseq = seqKmers[idx]
        if idx in revcomp_kmers:
            kseq = revcomp(kseq)
        filt[id] = kseq
        kmer_presence = []
        for seq_id in kmer_counts:
            if kmer_counts[seq_id][idx] == 1:
                kmer_presence.append(seq_id)
        kmer_presence = sorted(kmer_presence)
        fh.write("{}\t{}\t{}\t{}\t{}\n".format(id,s,e,kseq,",".join([str(x) for x in kmer_presence])))
        del(seqKmers[idx])
        id+=1
    return filt









