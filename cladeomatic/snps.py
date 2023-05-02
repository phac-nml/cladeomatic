import ray
from cladeomatic.utils.vcfhelper import vcfReader
from cladeomatic.utils import fisher_exact

@ray.remote
def snp_search(group_data, vcf_file, assigned_row, offset=0):
    """
    Accepts SNP data and tree to identify which SNPs correspond to a specific node on the tree.
    Following the search with :meth:`process_snps`, the data dictionary of snps is returned with
    group memberships and all other associated data for the snps is returned.

    Parameters
    ----------
    group_data : dict
        A dictionary of the snps and the clade/group membership
    ete_tree_object : ETE3 tree object
        The ETE3 object to traverse to associate snps to nodes
    vcf_file : str
        The file path to the VCF file or the TSV snp data file
    offset : int
        The amount to offset the rows.  Default of 0.

    Returns
    -------
    dict
        The dictionary of the :attr:`snp_data` data structure following the search through the available data
    Notes
    -----
    Refer to https://www.ray.io for more information about the Ray instances used in this module
    and http://etetoolkit.org for more thorough ETE3 documentation.

    """

    #read the vcf file
    vcf = vcfReader(vcf_file)
    data = vcf.process_row()
    samples = vcf.samples
    snps = {}
    count_snps = 0
    if data is not None:
        count_snps += 1
    #retrieve the dictionary of the sample ids, node ids and genotypes
    sample_map = group_data['sample_map']
    sample_id_lookup = {}
    for id in sample_map:
        sample_id_lookup[sample_map[id]['sample_id']] = id

    all_samples = set(sample_map.keys())
    #retrieve the dictionary of nodeIds and their clade/group members
    group_membership = group_data['membership']
    id = 0
    #loop through data to get the snp lists for the chromosomes in the data
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
            #assign the bases: group the ambiguous members and regular assignments
            ambig_members = set()
            if 'N' in assignments or '-' in assignments:
                if 'N' in assignments:
                    ambig_members = ambig_members | set(assignments['-']) | set(assignments['N'])
            for base in assignments:
                if not base in ['A', 'T', 'C', 'G']:
                    continue
                in_samples = set(assignments[base])
                is_ref = base == data['REF']
                #process the the snps - see method for documentation
                snps[chrom][pos][base] = process_snp(chrom, pos, base, all_samples, in_samples, ambig_members,
                                                     group_membership, is_ref)
            assigned_row += offset
        data = vcf.process_row()
        id += 1

    return snps


def snp_search_controller(group_data, vcf_file, n_threads=1):
    """
    A method to search for the lists of snps in the VCF file
    and group data to match them to create a complete snp dictionary.
    Please refer to the file examples/small_test/cladeomatic/cladeomatic-snps.info.txt
    for more information.

    Parameters
    ----------
    group_data : dict
        A dictionary of all the SNPs and their clade/group membership
    vcf_file : str
        The file path to the user VCF file containing the variants of interest
    n_threads : int
        The number of threads to use in the Ray method call

    Returns
    -------
    dict
        The data dictionary for the snps as identified from the input and grouped by chromosome, position and base

    Notes
    -----
    Refer to https://www.ray.io for more information about the Ray instances used in this module.

    """
    offsets = range(1, n_threads + 1)
    starts = range(0, n_threads)
    result_ids = []
    #loop through the number of threads to start ray instances
    #for the retrieval of snp value dictionaries
    for idx, value in enumerate(offsets):
        offset = value
        assigned_row = starts[idx]
        result_ids.append(
            snp_search.remote(group_data, vcf_file, assigned_row, offset))
    #retrieve the results from the snp_search method
    results = ray.get(result_ids)
    snps = {}
    #create the snps dictionary from the results
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
    """
    A method to process all the inputs to create a dictionary entry of SNP
    data.  Please refer to the file
    examples/small_test/cladeomatic/cladeomatic-snps.info.txt for more information.

    Parameters
    ---------
    chrom : str
        The chromosome to which the SNP belongs to
    pos : int
        The location of the SNP in the sequence
    base : str
        The nucleotide base of th SNP
    all_samples : set
        The set of identifiers for the samples
    snp_members : set
        The identifiers for the SNPs
    ambig_members : set
        The collections of ambiguously called nucleotides
    group_membership : dict
        The dictionary for the node ids and their clade/group memberships
    is_ref : bool
        True, if this the reference sequence nucleotide
    Returns
    _______
    dict
        A dictionary of the snp and associated information for further processing

    """
    num_members = len(snp_members)
    #create a set of all the samples without the ambiguous members
    all_samples = all_samples - ambig_members
    #create a set of snps without the ambiguous members
    snp_members = snp_members - ambig_members
    best_oddsr = 0
    best_p = 1
    best_clade_id = -1
    best_clade_num = 0
    is_canonical = False
    num_clade_members = 0

    #iterate through the group membership clades
    for clade_id in group_membership:
        if is_canonical:
            break
        #the members for the clade Id without the ambiguous members Ids
        clade_members = group_membership[clade_id] - ambig_members
        #find the members common to both the snp set and the clade set
        pos_pos = snp_members & clade_members
        #find the clade members NOT in the snp set
        neg_pos = clade_members - snp_members
        num_pos_pos = len(pos_pos)

        # Heuristic to skip poor quality comparisons
        if num_pos_pos == 0 or num_pos_pos / num_members < 0.5:
            continue
        #find the snps NOT in the clade memberships
        pos_neg = snp_members - clade_members
        #find the samples NOT in either the snps or the clade sets
        neg_neg = all_samples - snp_members - clade_members
        #create a table of the sets
        table = [[num_pos_pos, len(pos_neg | neg_neg)],
                 [len(neg_pos | pos_pos), len(neg_neg | pos_pos)]
                 ]
        #calculate the fisher's exact test for the above table
        oddsr, p = fisher_exact(table, alternative='greater')
        #determine the best clade id in the group membership with the best fisher score
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
