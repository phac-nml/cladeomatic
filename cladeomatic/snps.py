import ray
from cladeomatic.utils.vcfhelper import vcfReader
from cladeomatic.utils import fisher_exact

@ray.remote
def snp_search(group_data, vcf_file, assigned_row, offset=0):
    '''
    Accepts SNP data and tree to identify which SNPs correspond to a specific node on the tree
    :param ete_tree_obj: ETE3 tree object
    :param vcf_file: str path to vcf or tsv snp data
    :return: dict of snp_data data structure
    '''
    vcf = vcfReader(vcf_file)
    data = vcf.process_row()
    samples = vcf.samples
    snps = {}

    sample_map = group_data['sample_map']
    sample_id_lookup = {}
    for id in sample_map:
        sample_id_lookup[sample_map[id]['sample_id']] = id

    all_samples = set(sample_map.keys())
    group_membership = group_data['membership']
    id = 0

    while data is not None:
        if id == assigned_row:
            chrom = data['#CHROM']
            pos = int(data['POS']) - 1
            if not chrom in snps:
                snps[chrom] = {}
            snps[chrom][pos] = {}
            assignments = {}
            for sample_id in samples:
                base = data[sample_id]
                if not base in assignments:
                    assignments[base] = []
                assignments[base].append(sample_id_lookup[sample_id])

            ambig_members = set()
            if 'N' in assignments or '-' in assignments:
                if 'N' in assignments:
                    ambig_members = ambig_members | set(assignments['-']) | set(assignments['N'])
            for base in assignments:
                if not base in ['A', 'T', 'C', 'G']:
                    continue
                in_samples = set(assignments[base])
                is_ref = base == data['REF']
                snps[chrom][pos][base] = process_snp(chrom, pos, base, all_samples, in_samples, ambig_members,
                                                     group_membership, is_ref)
            assigned_row += offset
        data = vcf.process_row()
        id += 1

    return snps


def snp_search_controller(group_data, vcf_file, n_threads=1):
    '''

    Parameters
    ----------
    group_data
    vcf_file
    n_threads

    Returns
    -------

    '''
    offsets = range(1, n_threads + 1)
    starts = range(0, n_threads)
    result_ids = []
    for idx, value in enumerate(offsets):
        offset = value
        assigned_row = starts[idx]
        result_ids.append(
            snp_search.remote(group_data, vcf_file, assigned_row, offset))

    results = ray.get(result_ids)
    snps = {}
    for i in range(0, len(results)):
        result = results[i]
        for chrom in result:
            if not chrom in snps:
                snps[chrom] = {}
            for pos in result[chrom]:
                if not pos in snps[chrom]:
                    snps[chrom][pos] = {}
                for base in result[chrom][pos]:
                    snps[chrom][pos][base] = result[chrom][pos][base]
    return snps


def process_snp(chrom, pos, base, all_samples, snp_members, ambig_members, group_membership, is_ref):
    '''

    Parameters
    ----------
    chrom
    pos
    base
    all_samples
    snp_members
    ambig_members
    group_membership
    is_ref

    Returns
    -------

    '''
    num_members = len(snp_members)
    all_samples = all_samples - ambig_members
    snp_members = snp_members - ambig_members
    best_oddsr = 0
    best_p = 1
    best_clade_id = -1
    best_clade_num = 0
    is_canonical = False
    num_clade_members = 0

    for clade_id in group_membership:
        if is_canonical:
            break
        clade_members = group_membership[clade_id] - ambig_members
        pos_pos = snp_members & clade_members
        neg_pos = clade_members - snp_members
        num_pos_pos = len(pos_pos)

        # Heuristic to skip poor quality comparisons
        if num_pos_pos == 0 or num_pos_pos / num_members < 0.5:
            continue

        pos_neg = snp_members - clade_members
        neg_neg = all_samples - snp_members - clade_members

        table = [[num_pos_pos, len(pos_neg | neg_neg)],
                 [len(neg_pos | pos_pos), len(neg_neg | pos_pos)]
                 ]

        oddsr, p = fisher_exact(table, alternative='greater')



        if (oddsr > best_oddsr or p < best_p) and len(pos_neg) == 0:
            best_clade_id = clade_id
            best_clade_num = len(group_membership[clade_id])
            best_p = p
            best_oddsr = oddsr
            num_clade_members = len(clade_members)


        if len(pos_neg) == 0 and len(neg_pos) == 0:
            is_canonical = True
            best_clade_id = clade_id
            best_clade_num = len(group_membership[clade_id])
            best_p = p
            best_oddsr = oddsr
            num_clade_members = len(clade_members)


    return {'chrom': chrom, 'pos': pos, 'base': base,
            'clade_id': best_clade_id, 'is_canonical': is_canonical,'is_valid':True,
            'num_clade_members': num_clade_members, 'num_members': best_clade_num, 'is_ref': is_ref,
            'oddsr': best_oddsr, 'p_value': best_p}
