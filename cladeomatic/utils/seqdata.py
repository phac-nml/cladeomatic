import sys

from Bio import SeqIO, GenBank
import random, hashlib, copy, re
from cladeomatic.utils.vcfhelper import vcfReader


NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')

def revcomp(s):
    """Reverse complement nucleotide sequence

    Args:
        s (str): nucleotide sequence

    Returns:
        str: reverse complement of `s` nucleotide sequence
    """
    return s.translate(NT_SUB)[::-1]

def read_fasta_dict(fasta_file):
    """

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
    '''

    :param gbk_file: GenBank formatted sequence file
    :return: [dict] of sequences indexed by sequence id
    '''
    seqs = dict()
    with open(gbk_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            seqs[str(record.id)] = str(record.seq).upper()
    handle.close()
    return

def parse_reference_gbk(gbk_file):
    """
    :param gbk_file: Reference genbank format file with sequence annotations
    :return: dict of all of the reference features
    """
    sequences = {}
    with open(gbk_file) as handle:
        for record in GenBank.parse(handle):
            gb_accession = record.accession[0]
            gb_accession_version = gb_accession[1]
            genome_seq = repr(record.sequence).replace("\'",'')
            sequences[gb_accession] = {
                'accession':gb_accession,
                'version': gb_accession_version,
                'features': {'source': genome_seq}
            }
            features = record.features
            for feat in features:
                if feat.key == 'CDS' or feat.key == '5\'UTR' or feat.key == '3\'UTR':
                    if not feat.key in sequences[gb_accession]['features']:
                        sequences[gb_accession]['features'][feat.key] = []
                    qualifier = feat.qualifiers
                    positions = []
                    gene_name = ''
                    locus_tag = ''
                    aa = ''
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

                    for location in locations:
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
    '''
    :param string: string to comput MD5
    :return: md5 hash
    '''
    seq = str(string).encode()
    md5 = hashlib.md5()
    md5.update(seq)
    return md5.hexdigest()

def generate_non_gap_position_lookup(seq):
    """
    Creates a list of positions which correspond to the position of that base in a gapless sequence
    :param seq: string
    :return: list
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
    vcf = vcfReader(vcf_file)
    data = vcf.process_row()
    samples = vcf.samples
    valid_bases = ['A', 'T', 'C', 'G', '-', 'N']

    sample_variants = {}
    for sample in samples:
        sample_variants[sample] = {}

    if data is None:
        return {}

    while data is not None:
        chrom = data['#CHROM']
        pos = int(data['POS'])
        ref = data['REF']

        for sample_id in samples:
            base = data[sample_id]

            if base == '*':
                base = '-'
            if base not in valid_bases:
                base = 'N'
            is_ref = base == ref
            if is_ref:
                continue
            if not chrom in sample_variants[sample_id]:
                sample_variants[sample_id][chrom] = {}

            sample_variants[sample_id][chrom][pos] = base
        data = vcf.process_row()
    return sample_variants



def create_pseudoseqs_from_vcf(ref_seq,vcf_file, outfile):
    sample_variants = get_variants(vcf_file)
    fh = open(outfile,'w')
    for chr in ref_seq:
        fh.write(">{}~{}\n{}\n".format(chr, chr, ''.join(ref_seq[chr])))
    seqLens = {}
    chrom_id_map = {}
    id = 1
    for chrom in ref_seq:
        seqLens[chrom] = len(ref_seq[chrom])
        chrom_id_map[str(id)] = chrom
        id+=1


    for sample_id in sample_variants:
        if len(sample_variants[sample_id]) ==0:
            sample_variants[sample_id][list(chrom_id_map.keys())[0]] = {}
        for chrom in sample_variants[sample_id]:
            c = chrom
            if chrom not in ref_seq:
                if chrom in chrom_id_map:
                    c = chrom_id_map[chrom]
            seq = list(copy.deepcopy(ref_seq[c]))
            for pos in sample_variants[sample_id][chrom]:
                if pos > seqLens[c]:
                    print("Error variant position is outside sequence, check sequence for insertions which are not supported: {} seqlen {} pos".format(seqLens[c],pos))
                base = sample_variants[sample_id][chrom][pos]
                seq[pos-1] = base
            seq = ''.join(seq)
            fh.write(">{}~{}\n{}\n".format(sample_id, chrom, seq))

    fh.close()

def calc_homopolymers(seq):
    longest = 0
    for b in ['A', 'T', 'C', 'C']:
        matches = re.findall("{}+".format(b), seq)
        for m in matches:
            length = len(m)
            if length > longest:
                longest = length
    return longest