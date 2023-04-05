import sys

from Bio import SeqIO, GenBank
import random, hashlib, copy, re
from cladeomatic.utils.vcfhelper import vcfReader


NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')

def revcomp(s):
    """
    Reverse complement nucleotide sequence

    Args:
        s (str): nucleotide sequence

    Returns:
        str: reverse complement of `s` nucleotide sequence
    """
    return s.translate(NT_SUB)[::-1]

def read_fasta_dict(fasta_file):
    """
    Reads the fasta file from the passed file path and formats
    the input to a dictionary of sequences
    :param fasta_file: [str] Path to fasta file to read
    :return: [dict] of sequences indexed by sequence id
    """
    seqs = dict()
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs[str(record.id)] = str(record.seq).upper()
    handle.close()
    return seqs

def gb_to_fasta_dict(gbk_file):
    """
    Reads a GenBank formatted sequence file and creates a
    dictionary of sequences with the sequence id as keys
    :param gbk_file: String - the string path to a GenBank formatted
    sequence file
    :return: dictionary - a dictionary of sequences indexed by sequence id
    """
    seqs = dict()
    with open(gbk_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            seqs[str(record.id)] = str(record.seq).upper()
    handle.close()
    return

def parse_reference_gbk(gbk_file):
    """
    Method to parse the GenBank reference file, clean the strings, and
    return the reference features of interest.
    :param gbk_file: String - path to the reference genbank
    format file with sequence annotations
    :return: dictionary - a dictionary of all of the reference features
    """
    sequences = {}
    with open(gbk_file) as handle:
        for record in GenBank.parse(handle):
            gb_accession = record.locus
            gb_accession_version = 1
            if len(record.accession) == 2:
                gb_accession = record.accession[0]
                gb_accession_version = gb_accession[1]
            #clean the sequence
            genome_seq = repr(record.sequence).replace("\'",'')
            sequences[gb_accession] = {
                'accession':gb_accession,
                'version': gb_accession_version,
                'features': {'source': genome_seq}
            }
            features = record.features
            #retrieve the features if present
            for feat in features:
                if feat.key == 'CDS' or feat.key == '5\'UTR' or feat.key == '3\'UTR':
                    if not feat.key in sequences[gb_accession]['features']:
                        sequences[gb_accession]['features'][feat.key] = []
                    qualifier = feat.qualifiers
                    positions = []
                    gene_name = ''
                    locus_tag = ''
                    aa = ''
                    #more string cleaning
                    for name in qualifier:
                        if name.key == '/gene=':
                            gene_name = name.value.replace("\"", '').strip()
                        if name.key == '/translation=':
                            aa = name.value.replace("\"", '').strip()
                        if name.key == '/locus_tag=':
                            gene_name = name.value.replace("\"", '').strip()
                            locus_tag = gene_name
                    if locus_tag != '':
                        gene_name = locus_tag
                    locations = feat.location.strip().replace("join(", '').replace(')', '').split(',')
                    seq = []
                    #retreive the locations of the features
                    for location in locations:
                        #more string cleaning
                        location = location.replace('<','').replace('>','')
                        if not 'complement' in location:
                            location = location.split('.')
                            start = int(location[0]) - 1
                            end = int(location[2])
                            seq.append(genome_seq[start:end].replace("\'", ''))
                            positions.append([start, end])
                        else:
                            location = location.replace('complement(','').replace(')','').split('.')
                            start = int(location[0]) - 1
                            end = int(location[2])
                            seq.append(revcomp(genome_seq[start:end].replace("\'", '')))
                            positions.append([start, end])

                    seq = ''.join(seq)
                    sequences[gb_accession]['features'][feat.key].append(
                        {'gene_name': gene_name, 'dna_seq': seq, 'aa_seq': aa, 'positions': positions,'gene_len':len(seq)})

    return sequences

def calc_md5(string):
    """
    Method to encode the MD5 hash for the input string
    :param string: String to compute MD5
    :return: hash - the md5 hash
    """
    seq = str(string).encode()
    md5 = hashlib.md5()
    md5.update(seq)
    return md5.hexdigest()

def generate_non_gap_position_lookup(seq):
    """
    Creates a list of positions which correspond to the position
    of that base in a gapless sequence
    :param seq: string - the sequence to process
    :return: list - an int list of the positions  of the gaps
    """
    length = len(seq)
    num_gaps = 0
    lookup = []
    for i in range(0, length):
        base = seq[i]
        if base == '-':
            num_gaps += 1
            lookup.append(-1)
        else:
            lookup.append(i - num_gaps)
    return lookup

def create_aln_pos_from_unalign_pos_lookup(aln_seq):
    """
    This method creates a list of integers for the positions
    of the bases in the unaligned sequence derived from the
    aligned sequence passed to the method
    :param aln_seq: String - the alignment sequence
    :return: list - a list of integers for the positions found
    """
    unalign_seq = aln_seq.replace('-', '')
    aln_len = len(aln_seq)
    unaln_len = len(unalign_seq)
    lookup = [-1] * unaln_len
    pos = 0
    for i in range(0, unaln_len):
        for k in range(pos, aln_len):
            if unalign_seq[i] == aln_seq[k]:
                lookup[i] = k
                pos = k + 1
                break
    return lookup

def get_variants(vcf_file):
    """
    This method reads the incoming VCF file and returns
    a dictionary of the sample variant bases and their locations
    both in the genome and the tree
    :param vcf_file: String - path to the VCF file
    :return: dictionary - a dictionary of the sample variants with the
    node id as key
    """
    #read the file
    vcf = vcfReader(vcf_file)
    data = vcf.process_row()
    samples = vcf.samples
    #valid bases include missing data - still valid for our purposes
    valid_bases = ['A', 'T', 'C', 'G', '-', 'N']

    sample_variants = {}
    for sample in samples:
        sample_variants[sample] = {}

    if data is None:
        return {}
    #process the VCF data
    while data is not None:
        chrom = data['#CHROM']
        pos = int(data['POS']) -1
        ref = data['REF']
        #link the samples id and the chromosomal location to their variant base
        for sample_id in samples:
            base = data[sample_id]

            if base == '*':
                base = '-'
            if base not in valid_bases:
                base = 'N'
            is_ref = base == ref
            #if the base is the same as the reference base, skip
            if is_ref:
                continue
            if not chrom in sample_variants[sample_id]:
                sample_variants[sample_id][chrom] = {}

            sample_variants[sample_id][chrom][pos] = base
        data = vcf.process_row()
    return sample_variants

def create_pseudoseqs_from_vcf(ref_id,ref_seq,vcf_file, outfile):
    """
    This method creates the pseudo full sequences of the variants found
    in the VCF file.
    :param ref_id: String - the reference sequence identifier
    :param ref_seq: String - the reference sequence to alter
    :param vcf_file: String - the path to the VCF file
    :param outfile: String - the path to the pseudo variant outfile
    """
    #retrieve the variant bases, genome and tree locations
    sample_variants = get_variants(vcf_file)
    fh = open(outfile,'w')
    fh.write(">{}\n{}\n".format(ref_id, ref_seq))
    ref_len = len(ref_seq)
    #loop through the sample variants to replace the reference base
    #with the variant base at the variant base position
    #to create the pseudo sequence
    for sample_id in sample_variants:
        if sample_id == ref_id:
            continue
        if len(sample_variants[sample_id]) ==0:
            sample_variants[sample_id] = {}
        #deep copy the reference sequence
        seq = list(copy.deepcopy(ref_seq))
        #perform the replacement for the variant
        for chrom in sample_variants[sample_id]:
            for pos in sample_variants[sample_id][chrom]:
                if pos >= ref_len:
                    print("Error variant position is outside sequence, check sequence for insertions which are not supported: {} seqlen {} pos".format(ref_len,pos))
                    sys.exit()
                base = sample_variants[sample_id][chrom][pos]
                seq[pos] = base
            seq = ''.join(seq)
            fh.write(">{}\n{}\n".format(sample_id,seq))
    fh.close()

def calc_homopolymers(seq):
    """
    The method calculates the longest homopolymer (the sequence
    of consecutive identical bases) in the sequence or sequence
    fragment passed
    :param seq: String - the sequences or sequence fragment to find the
    longest homopolymer
    :return: int - the longest homopolymer length
    """
    longest = 0
    for b in ['A', 'T', 'C', 'C']:
        matches = re.findall("{}+".format(b), seq)
        for m in matches:
            length = len(m)
            if length > longest:
                longest = length
    return longest