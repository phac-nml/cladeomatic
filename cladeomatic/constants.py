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


GENOTYPE_REPORT_HEADER = [
    'sample_id',
    'scheme',
    'analysis_date',
    'total_scheme_features',
    'unique_scheme_positions',
    'num_detected_features',
    'num_detected_positions',
    'predicted_genotype',
    'mutations_supporting_prediction',
    'confidence',
    'qc_status',
    'qc_messages',
    'detected_alt_mutations'
]

NODE_INFO_HEADER = [
        'clade_id',
        'pos',
        'base',
        'min_dist',
        'max_dist',
        'ave_dist',
        'num_members',
        'is_valid',
        'is_selected',
        'total_within_clade_dist',
        'num_comparisons_within_clade_dist',
        'ave_within_clade_dist',
        'closest_clade_id',
        'closest_clade_dist',
        'closest_sample_id',
        'closest_sample_dist',
        'clade_sample_id',
        'spearmanr',
        'spearmanr_pvalue',
        'pearsonr',
        'pearsonr_pvalue',
        'is_temporal_signal_present',
        'field_name',
        'fisher_oddsr',
        'fisher_p'
]


CONFIDENCE = ['strong','moderate','weak']

QC_STATUS = ['PASS','FAIL','WARN']

MIN_FILE_SIZE = 32



IUPAC_LOOK_UP = {
    'A': 'A',
    'T': 'T',
    'C': 'C',
    'G': 'G',
    'AT': 'W',
    'TA': 'W',
    'AC': 'M',
    'CA': 'M',
    'AG': 'R',
    'GA': 'R',
    'CT': 'Y',
    'TC': 'Y',
    'GT': 'K',
    'TG': 'K',
    'CG': 'S',
    'GC': 'S',
    'CGT': 'B',
    'GTC': 'B',
    'TCG': 'B',
    'CTG':'B',
    'TGC': 'B',
    'AGT': 'D',
    'GTA': 'D',
    'TAG': 'D',
    'GAT': 'D',
    'ATG': 'D',
    'ACT': 'H',
    'CTA': 'H',
    'TCA': 'H',
    'TAC': 'H',
    'CAT': 'H',
    'ACG': 'V',
    'CGA': 'V',
    'GCA': 'V',
    'CAG': 'V',
    'GAC': 'V',
    'ACGT': 'N',

}