LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

SCHEME_HEADER = [
    'key',
    'mutation_key',
    'mutation_type',
    'dna_name',
    'gene_name',
    'gene_start',
    'gene_end',
    'cds_start',
    'cds_end',
    'aa_name',
    'aa_start',
    'aa_end',
    'is_silent',
    'is_cds',
    'is_frame_shift',
    'variant_start',
    'variant_end',
    'kmer_start',
    'kmer_end',
    'target_variant',
    'target_variant_len',
    'state',
    'ref_state',
    'alt_state',
    'kseq',
    'klen',
    'homopolymer_len',
    'kmer_entropy',
    'positive_genotypes',
    'partial_genotypes',
    'is_ambig_ok',
    'is_kmer_found',
    'is_kmer_length_ok',
    'is_kmer_unique',
    'is_valid'
]

CLADE_DEFINING_SNPS_HEADER = [
    'genotype',
    'pos',
    'base',
    'dna_name',
    'is_cds',
    'gene_name',
    'aa_name'
    'is_silent',
]




IUPAC_LOOK_UP = {
        'A': 'A',
        'T': 'T',
        'C': 'C',
        'G': 'G',
        'AT': 'W',
        'AC': 'M',
        'AG': 'R',
        'CT': 'Y',
        'GT': 'K',
        'CG': 'S',
        'CGT': 'B',
        'AGT': 'D',
        'ACT': 'H',
        'ACG': 'V',
        'ACGT': 'N'
    }