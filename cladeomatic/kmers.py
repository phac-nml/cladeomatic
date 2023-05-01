import sys
import time

from Bio import SeqIO

from cladeomatic.utils.kmerSearch import SeqSearchController
from cladeomatic.constants import IUPAC_LOOK_UP
from deprecated import deprecated

class kmer_worker:
    """  The kmer_worker class.

    The kmer_worker class contains several class variables for
    use in the creation of the kmer lists required for downstream processes in the
    creation of the schemes and data analysis files.
    """
    # input
    ref_sequence = ''
    target_positions = []
    msa_fasta_file = None
    kmer_len = 0
    num_threads = 1
    result_dir = None
    prefix = 'cladeomatic'
    genotype_membership = {}
    max_ambig = 0

    # Derived
    ref_len = 0
    min_geno_perc = 1
    valid_bases = ['A', 'T', 'C', 'G']
    variant_positions = []
    msa_base_counts = {}
    opt_kmer_start_positions = {}
    extracted_kmers = {}
    num_extracted_kmers = 0
    kmer_search_files = []
    genotype_snp_associations = {}
    kmer_scheme_data = {}
    rule_set = {}
    invalid_kmer_indexes = set()
    int_base_kmer_lookup = {}
    positions_missing_kmer = {}
    biohansel_kmers = {}
    genotype_snp_rules = {}

    def __init__(self, ref_sequence, msa_file, result_dir, prefix, klen, genotype_map, genotype_snp_rules,max_ambig=0, min_perc=1,
                 target_positions=[], num_threads=1):
        """
        The instantiation of the :class:`kmer_worker` class.

        Parameters
        ----------
        ref_sequence : str
            The reference sequence for kmer selection
        msa_file : str
            The file path to the fasta file with the snps substitutions (multiple sequence alignment)
        result_dir : str
            The file path to the results directory
        prefix : str
            The prefix for the clade-o-matic produced output files
        klen : int
            The length of the kmers
        genotype_map : dict
            The map of the sample identifiers and genotypes
        genotype_snp_rules : dict
            The dictionary groupings of the positions, base variants, positive and partial genotypes for the snps, that make up the processing rules
        max_ambig : int
            The maximum number of ambiguous bases that can be contained in a kmer sequence.  Default for this class is 0.
        min_perc : float
            The minimum percentage of clade members to be positive for a kmer to be considered valid. Default for this class is 1.
        target_positions : list
            The integer list of the target snp positions.  Default for this class is an empty list.
        num_threads : int
            The number of threads for the Ray instance.  Default for this class is 1.

        Notes
        -----
        Refer to https://www.ray.io for more information about the Ray instances used in this module
        """
        self.ref_sequence = ref_sequence
        self.ref_len = len(ref_sequence)
        self.max_ambig = max_ambig
        self.msa_fasta_file = msa_file
        self.target_positions = target_positions
        self.num_threads = num_threads
        self.kmer_len = klen
        self.min_geno_perc = min_perc
        self.result_dir = result_dir
        self.prefix = prefix
        self.genotype_membership = genotype_map
        self.genotype_snp_rules = genotype_snp_rules
        self.workflow()
        return

    def workflow(self):
        """
        The workflow method to call all the methods within this class to create the kmer lists
        used for schema analysis and processing.
        """
        self.init_msa_base_counts()
        self.pop_msa_base_counts()
        self.find_variant_positions()
        if len(self.target_positions) == 0:
            self.target_positions = self.variant_positions
        self.get_genotype_snp_states()
        self.get_optimal_kmer_position()
        self.extract_kmers()
        self.populate_int_base_kmer_lookup()
        self.perform_kmer_search()
        self.process_kmer_results()
        self.init_kmer_scheme_data()
        self.populate_kmer_scheme_data()
        self.confirm_kmer_specificity()


        # Reinitialize with invalid kmers removed
        self.init_kmer_scheme_data()
        self.populate_kmer_scheme_data()
        self.remove_empty_base_states()
        self.construct_ruleset()
        self.remove_redundant_kmers()
        self.find_invalid_kmers()
        self.remove_invalid_kmers_from_scheme()
        self.positions_missing_kmer = self.get_pos_without_kmer()
        self.refine_rules()
        self.biohansel_kmers = self.create_biohansel_kmers()


    def init_msa_base_counts(self):
        """
        Method to initialize the msa (multiple sequence alignment) base counts dictionary :attr:`msa_base_counts`
        with the number of each base substitution for the length of the reference sequence set to
        zero.
        """
        for i in range(0, self.ref_len):
            self.msa_base_counts[i] = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, '-': 0}

    def pop_msa_base_counts(self):
        """
        Method to count or populate the :attr:`msa_base_counts` dictionary initialized in
        :meth:`init_msa_base_counts` for the bases in the fasta file passed.
        For the sequence in the fasta, the snp position is referenced and the counter
        is updated for each base substitution that occurs in that position.
        """
        with open(self.msa_fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper()
                for pos in self.msa_base_counts:
                    base = seq[pos]
                    if not base in self.msa_base_counts[pos]:
                        base = 'N'
                    self.msa_base_counts[pos][base] += 1
            handle.close()
    def find_variant_positions(self):
        """
        Method to find the sequence variants and their positions in the
        :attr:`msa_base_counts` dictionary and assign the positions of the bases
        to the variant_positions dictionary.
        """
        for pos in self.msa_base_counts:
            count_bases = 0
            #ensure the base is valid and if there is more than one, add it
            #as a variant
            for base in self.valid_bases:
                if self.msa_base_counts[pos][base] > 0:
                    count_bases += 1
            if count_bases > 1:
                self.variant_positions.append(pos)

    def get_optimal_kmer_position(self):
        """
        A method to determine the best kmer start positions for the
        variant/snps positions found in the find_variant_positions method.
        Set these kmer start positions in a list for further processing
        """
        msa_len = self.ref_len
        klen = self.kmer_len
        local_kmer_start_positions = set()
        #loop through the variant positions
        for pos in self.target_positions:
            start = pos - klen + 1
            if start < 0:
                start = 0
            end = pos + 1
            if end > msa_len:
                end = msa_len
            interval = range(start, end)
            best_count_variable_sites = klen
            best_count_missing_sites = klen
            best_index = start
            #determine the best start position for the kmer of the
            #variant based on the reference sequence
            for i in interval:
                variable_sites = 0
                num_missing = 0

                for k in range(i, i + klen + 1):
                    if k >= end:
                        break
                    num_bases = 0
                    # determine if there are missing bases in the kmer
                    if self.msa_base_counts[k]['-'] > 0:
                        num_missing += 1

                    for base in self.valid_bases:
                        if self.msa_base_counts[k][base] > 0:
                            num_bases += 1
                    if num_bases > 1:
                        variable_sites += 1
                #the logic for choosing the best kmers
                #if there is a snp present and the selected kmer has no missing bases,
                #it is the best representative kmer
                if variable_sites <= best_count_variable_sites:
                    if num_missing < best_count_missing_sites:
                        best_index = i
                        best_count_variable_sites = variable_sites
                        best_count_missing_sites = num_missing

                local_kmer_start_positions.add(best_index)
        self.opt_kmer_start_positions = sorted(list(local_kmer_start_positions))

    def extract_kmers(self):
        """
        This method creates the dictionary of all the possible kmers
        for the snp variants found in the sequence alignments.  These kmers are
        found through the looping of the snp positions in the sequence and 'frame-shifting'
        the kmer sequence through all possibilities for the snp, as long as there are no
        ambiguous bases nor gaps in the kmer sequence and the kmer is the specified length.
        """
        msa_len = self.ref_len
        klen = self.kmer_len
        selected_kmers = {}
        canonical_pos = set(self.target_positions)
        #read the fasta
        with open(self.msa_fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper()
                #loop through the snp locations
                for pos in self.target_positions:
                    base = seq[pos]
                    if base not in self.valid_bases:
                        continue
                    #get the starting position for the kmer
                    s = pos - klen + 1
                    if s < 0:
                        s = 0
                    start_range = range(s, pos)
                    #create a range of potential kmers for the snp based on the start position
                    for start_pos in start_range:
                        start = start_pos

                        end = start + klen
                        #extract the kmer from the sequence
                        kseq = seq[start:end + 1]
                        min_ambig = kseq.count('N')
                        no_gap_len = len(kseq.replace("-", ""))
                        # if the length of the kmer is longer than the specified kmer length
                        # AND the end position is larger than the position of the snp
                        # AND the number of ambiguous bases is less than or equal to the
                        # set minimum ambiguous level
                        # THEN shift the end less 1 base and remove ambiguous bases if present
                        while no_gap_len > klen and end > pos and kseq.count('N') <= min_ambig:
                            end -= 1
                            kseq = seq[start:end + 1]
                            no_gap_len = len(kseq.replace("-", ""))
                            if kseq.count('N') <= min_ambig:
                                min_ambig = kseq.count('N')
                        # if the length of the kmer is shorter than the specified kmer length
                        # AND the end position of the kmer is shorter than the sequence length
                        # AND the number of ambiguous bases is less than or equal to the
                        # set minimum ambiguous level
                        # Then add a base to the kmer and remove ambiguous bases if present
                        while no_gap_len < klen and end < msa_len and kseq.count('N') <= min_ambig:
                            end += 1
                            kseq = seq[start:end + 1]
                            no_gap_len = len(kseq.replace("-", ""))
                            if kseq.count('N') <= min_ambig:
                                min_ambig = kseq.count('N')
                        # if the length of the kmer is shorter than the specified kmer length
                        # AND the start position of the kmer is negative
                        # AND the number of ambiguous bases is less than or equal to the
                        # set minimum ambiguous level
                        # Then subtract a base to the kmer and remove ambiguous bases if present
                        while no_gap_len < klen and start > 0 and kseq.count('N') <= min_ambig:
                            start -= 1
                            kseq = seq[start:end + 1]
                            no_gap_len = len(kseq.replace("-", ""))
                            if kseq.count('N') <= min_ambig:
                                min_ambig = kseq.count('N')
                        #remove any ambiguous bases
                        kseq = seq[start:end + 1].replace('-', '')
                        if kseq.count('N') > self.max_ambig or len(kseq) != klen:
                            continue
                        #create the kmer dictionary
                        tpos = set(range(start, end)) & canonical_pos
                        tpos_bases = {}
                        for p in tpos:
                            tpos_bases[p] = seq[p]
                        selected_kmers[kseq] = {'aln_start': start, 'aln_end': end, 'is_valid': True,
                                                'target_positions': tpos_bases, 'genotype_counts': {}}

            handle.close()
        self.extracted_kmers = selected_kmers
        self.num_extracted_kmers = len(selected_kmers)

    def perform_kmer_search(self):
        """
        This method takes in the list of all possible kmers for the sequences
        provided from the method :meth:`extract_kmers` and writes them to a file for later
        processing using the :meth:`cladeomatic.utils.kmerSearch.SeqSearchController` method.
        """
        seqkmers = {}
        index = 0
        for kmer in self.extracted_kmers:
            seqkmers[index] = kmer
            index += 1
        self.kmer_search_files = SeqSearchController(seqkmers, self.msa_fasta_file, self.result_dir, self.prefix,
                                                     self.num_threads)

    def process_kmer_results(self):
        """
        This method process the all the kmers within the search file clade-o-matic created
        in the :meth:`perform_kmer_search` method.  This method filters and find the list
        of valid kmers for the sequence variants by ensuring the kmer sequence id exists
        in the genotype map and there are no duplicate kmers in the search file
        """
        invalid_kmers = set()
        genotype_mapping = self.genotype_membership
        #ensure the kmer sequence id exists in the genotype map
        for filename in self.kmer_search_files:
            fh = open(filename, 'r')
            for line in fh:
                line = line.rstrip().split("\t")
                seq_id = str(line[0]).split("~")[0]
                if not seq_id in genotype_mapping:
                    continue
                genoytpe = genotype_mapping[seq_id]
                kIndex = int(line[1])
                count = int(line[2])

                if count == 0:
                    continue
                #if there are more than 1 kmer in the search file,
                #add it to the invalid list
                if count > 1:
                    invalid_kmers.add(kIndex)
                if not genoytpe in self.int_base_kmer_lookup[kIndex]['genotype_counts']:
                    self.int_base_kmer_lookup[kIndex]['genotype_counts'][genoytpe] = 0
                self.int_base_kmer_lookup[kIndex]['genotype_counts'][genoytpe] += 1
            self.flag_invalid_kmers(invalid_kmers)

    def flag_invalid_kmers(self, invalid_kmers):
        """
        This method takes the set of invalid kmer indexes and flags the kmers as invalid in
        the extracted kmer dictionary.

        Parameters
        ----------
        invalid_kmers : set
            The set of integers for the indexes of the invalid kmers that require flagging in the extracted kmes dictionary
        """
        kIndex = 0
        for kmer in self.extracted_kmers:
            if kIndex in invalid_kmers:
                self.extracted_kmers[kmer]['is_valid'] = False
            kIndex += 1

    def init_kmer_scheme_data(self):
        """
        A method to initialize the kmer scheme through the list of
        target positions.  Adds blank entries to the kmer scheme
        dictionary.
        """
        kmer_data = {}
        for pos in self.target_positions:
            kmer_data[pos] = {
                'A': [], 'T': [], 'C': [], 'G': [],
            }
        self.kmer_scheme_data = kmer_data

    def populate_kmer_scheme_data(self):
        """
        A method to populate the kmer scheme dictionary with the
        processed kmers of the scheme data and extracted kmers dictionaries.
        This method removes the invalid kmers, removes the kmers with invalid bases,
        and the kmers with invalid positions in comparison to the scheme data.
        """
        scheme_data = self.kmer_scheme_data
        kmer_data = self.extracted_kmers
        index = 0
        invalid_kmers = set()
        for kmer in kmer_data:
            is_valid = kmer_data[kmer]['is_valid']
            ovl_pos = kmer_data[kmer]['target_positions']
            #if the kmer is not valid, skip it and increase the index
            if not is_valid:
                index += 1
                continue
            for pos in ovl_pos:
                base = ovl_pos[pos]
                if base not in self.valid_bases:
                    invalid_kmers.add(index)
                    continue
                if not pos in scheme_data:
                    continue
                scheme_data[pos][base].append(index)
            index += 1

    def get_pos_without_kmer(self):
        """
        This method compiles a dictionary of missing kmer indexes as identifiers
        and the snp base they flank.

        Returns
        -------
        dict
            A dictionary of the missing kmers and the snp bases these kmers flank.
        """
        missing = {}
        for pos in self.kmer_scheme_data:
            base_counts = self.msa_base_counts[pos]
            valid_bases = {}
            for b in self.valid_bases:
                if base_counts[b] > 0:
                    valid_bases[b] = base_counts[b]
                    if b not in self.kmer_scheme_data[pos] or len(self.kmer_scheme_data[pos][b]) == 0:
                        if not pos in missing:
                            missing[pos] = []
                        missing[pos].append(b)

        return missing

    def construct_ruleset(self):
        """
        A method to construct and populate the kmer rule set dictionary
        of the positive genotypes (the kmers that match or exceed the minimum percentage
        of clade members to be positive for a kmer to be valid) and the partial
        genotypes which are less than the minimum percentage of clade members.
        """
        genotype_map = self.genotype_membership
        #Minimum percentage of clade members to be positive for a kmer to be valid
        min_perc = self.min_geno_perc
        genotype_counts = {}
        for sample_id in genotype_map:
            genotype = genotype_map[sample_id]
            if not genotype in genotype_counts:
                genotype_counts[genotype] = 0
            genotype_counts[genotype] += 1

        kmer_info = self.extracted_kmers
        kmer_rules = {}
        index = 0
        for kmer in kmer_info:
            kmer_rules[index] = {'positive_genotypes': [], 'partial_genotypes': []}

            genotype_data = kmer_info[kmer]['genotype_counts']
            for genotype in genotype_data:
                if not genotype in genotype_counts:
                    continue
                total = genotype_counts[genotype]
                g = genotype_data[genotype]
                perc = g / total
                if perc >= min_perc:
                    kmer_rules[index]['positive_genotypes'].append(genotype)
                elif g > 0:
                    kmer_rules[index]['partial_genotypes'].append(genotype)
            index += 1


        self.rule_set = kmer_rules

    def refine_rules(self):
        """
        A method to filter and refine the kmer selection rules.  Missing data can cause
        kmers to all be partial when there isn't an alternative kmer available for a genotype.
        This filters the rules to assign kmers to be positive for a genotype when there
        is only one kmer base state present for it.
        """
        kmer_rules = self.rule_set
        #get all defined genotypes in the snp scheme
        genotypes_defined = set()
        for pos in self.genotype_snp_rules:
            for base in self.genotype_snp_rules[pos]:
                genotypes_defined = genotypes_defined | set(self.genotype_snp_rules[pos][base]['positive_genotypes']) \
                                    | set(self.genotype_snp_rules[pos][base]['partial_genotypes'])

        for pos in self.kmer_scheme_data:
            if pos == 0:
                continue
            positive_genos = {'A': set(), 'T': set(), 'C': set(), 'G': set()}
            partial_genos = {'A': set(), 'T': set(), 'C': set(), 'G': set()}
            target_genotypes = set()
            genos_with_positive_kmer = set()
            snp_rule = self.genotype_snp_rules[pos]
            for base in self.kmer_scheme_data[pos]:
                #create the sets for the positive and partial genotypes based on inclusion rules
                for kIndex in self.kmer_scheme_data[pos][base]:
                    positive_genos[base] = positive_genos[base] | set(kmer_rules[kIndex]['positive_genotypes'])
                    partial_genos[base] = partial_genos[base] | set(kmer_rules[kIndex]['partial_genotypes'])
                    target_genotypes = target_genotypes | positive_genos[base] |  partial_genos[base]
                    genos_with_positive_kmer = genos_with_positive_kmer | positive_genos[base]

                positive_genos[base] = set(positive_genos[base])
                positive_genos[base] = set(positive_genos[base])

            genos_missing_rule = genotypes_defined - target_genotypes
            for genotype in genos_missing_rule:
                for base in snp_rule:
                    if genotype in snp_rule[base]:
                        break
                if base in self.kmer_scheme_data[pos]:
                    if len(self.kmer_scheme_data[pos][base]) == 1:
                        genos_with_positive_kmer = genos_with_positive_kmer | set(genotype)
                        for kIndex in self.kmer_scheme_data[pos][base]:
                            kmer_rules[kIndex]['positive_genotypes'].append(genotype)
                    else:
                        for kIndex in self.kmer_scheme_data[pos][base]:
                            kmer_rules[kIndex]['partial_genotypes'].append(genotype)


            genotypes_to_check = target_genotypes - genos_with_positive_kmer
            for genotype in genotypes_to_check:
                bases_present = []
                for base in partial_genos:
                    if genotype in partial_genos[base]:
                        bases_present.append(base)
                if len(bases_present) > 1:
                    continue
                for base in bases_present:
                    affected_kmers = []

                    for kIndex in self.kmer_scheme_data[pos][base]:
                        if genotype in kmer_rules[kIndex]['partial_genotypes']:
                            affected_kmers.append(kIndex)

                    for kIndex in affected_kmers:
                            r = set(kmer_rules[kIndex]['partial_genotypes'])
                            r = r - set(genotype)
                            kmer_rules[kIndex]['partial_genotypes'] = list(r)
                            kmer_rules[kIndex]['positive_genotypes'].append(genotype)

        for kIndex in kmer_rules:
            kmer_rules[kIndex]['positive_genotypes'] = sorted(list(set(kmer_rules[kIndex]['positive_genotypes'])))
            kmer_rules[kIndex]['partial_genotypes'] = sorted(list(set(kmer_rules[kIndex]['partial_genotypes'])))

        self.rule_set = kmer_rules

    def get_genotype_snp_states(self):
        """
        A method to find the valid SNPs in the genotype membership dictionary
        previously constructed.
        """
        with open(self.msa_fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                id = str(record.id).split("~")[0]
                if id not in self.genotype_membership:
                    continue
                seq = str(record.seq).upper()
                #get the genotype id
                genotype = self.genotype_membership[id]
                #initialize the dictionary
                if not genotype in self.genotype_snp_associations:
                    self.genotype_snp_associations[genotype] = {}
                for pos in self.target_positions:
                    if not pos in self.genotype_snp_associations[genotype]:
                        self.genotype_snp_associations[genotype][pos] = set()
                    base = seq[pos]
                    #if the base is not valid for the position chosen, do not add to dictionary
                    if not base in self.valid_bases:
                        continue
                    #add the base change for the snp associated with the genotype id and sequence position
                    self.genotype_snp_associations[genotype][pos].add(base)
            handle.close()

    def confirm_kmer_specificity(self):
        """
        A method to further filter invalid kmers from the extracted kmer dictionary.
        This method ensures the target positions have the valid bases and exist in the
        kmer scheme data dictionary.  Kmers that do not adhere to these rules are flagged
        as invalid.
        """
        index = 0
        kmer_info = self.extracted_kmers
        invalid_kmers = set()
        for kmer in kmer_info:
            kmer_genotype_data = kmer_info[kmer]['genotype_counts']
            ovl_positions = kmer_info[kmer]['target_positions']
            for pos in ovl_positions:
                for genotype in kmer_genotype_data:
                    allowed_bases = self.genotype_snp_associations[genotype][pos]
                    is_found = False
                    for base in self.kmer_scheme_data[pos]:
                        if index in self.kmer_scheme_data[pos][base]:
                            is_found = True
                            break
                    if not is_found:
                        continue
                    if base not in allowed_bases:
                        invalid_kmers.add(index)
            index += 1
        self.flag_invalid_kmers(invalid_kmers)

    def find_invalid_kmers(self):
        """
        A method to loop through the extracted kmer dictionary,
        retrieve the flagged invalid kmers and add their indexes
        to a set of just invalid kmers.  This is a helper method used by other
        methods for kmer filtering.
        """
        index = 0
        for kmer in self.extracted_kmers:
            if self.extracted_kmers[kmer]['is_valid'] == False:
                self.invalid_kmer_indexes.add(index)
            index += 1

    def populate_int_base_kmer_lookup(self):
        """
        A method to instantiate the base kmer lookup dictionary
        with a sequential index and the sequence of the kmer.
        """
        index = 0
        for kmer in self.extracted_kmers:
            self.int_base_kmer_lookup[index] = self.extracted_kmers[kmer]
            index += 1

    @deprecated()
    def find_ovl_kmers(self):
        variant_positions = self.target_positions
        selected_kmers = self.kmer_scheme_data
        num_var_pos = len(variant_positions)
        shared_kmers = set()
        pot_pos_ovl = {}
        for i in range(0, len(variant_positions)):
            var_pos_1 = variant_positions[i]
            var_kmer_indices_1 = set()
            for base in selected_kmers[var_pos_1]:
                var_kmer_indices_1 = var_kmer_indices_1 | set(selected_kmers[var_pos_1][base])

            for k in range(i + 1, num_var_pos):
                var_pos_2 = variant_postions[k]
                if var_pos_2 - var_pos_1 > self.kmer_len:
                    break
                var_kmer_indices_2 = set()
                for base in selected_kmers[var_pos_2]:
                    var_kmer_indices_2 = var_kmer_indices_2 | set(selected_kmers[var_pos_2][base])

                shared_kmers = shared_kmers | (var_kmer_indices_1 & var_kmer_indices_2)
                if not var_pos_1 in pot_pos_ovl:
                    pot_pos_ovl[var_pos_1] = []
                pot_pos_ovl[var_pos_1].append(var_pos_2)
        return pot_pos_ovl

    def remove_invalid_kmers_from_scheme(self):
        """
        A method to remove the invalid kmers from the scheme data
        dictionary.
        """
        for index in self.kmer_scheme_data:
            for base in self.kmer_scheme_data[index]:
                #use the set data structure to remove the invalid
                #kmers from the kmer scheme dictionary
                ovl = set(self.kmer_scheme_data[index][base]) - self.invalid_kmer_indexes
                self.kmer_scheme_data[index][base] = ovl

    def remove_empty_base_states(self):
        """
        A clean up method to remove the records with empty or invalid bases from the kmer
        scheme data dictionary.
        """
        for index in self.kmer_scheme_data:
            for base in self.valid_bases:
                if base not in self.kmer_scheme_data[index]:
                    continue
                if len(self.kmer_scheme_data[index][base]) == 0:
                    del (self.kmer_scheme_data[index][base])

    def remove_redundant_kmers(self):
        """
        A method to remove all but the best kmer sequences from the kmer scheme
        data dictionary.  This is accomplished by revisiting the kmer rules sets
        to find the best representative kmer and finally updates the kmer data scheme
        dictionary with the results of the filtering.
        """
        kmer_info = self.int_base_kmer_lookup
        kmer_rules = self.rule_set
        bases = self.valid_bases
        num_bases = len(bases)
        for pos in self.kmer_scheme_data:
            if pos == 0:
                continue
            positive_genos = {'A': set(), 'T': set(), 'C': set(), 'G': set()}
            partial_genos = {'A': set(), 'T': set(), 'C': set(), 'G': set()}
            for base in self.kmer_scheme_data[pos]:
                for kIndex in self.kmer_scheme_data[pos][base]:
                    positive_genos[base] = positive_genos[base] | set(kmer_rules[kIndex]['positive_genotypes'])
                    partial_genos[base] = partial_genos[base] | set(kmer_rules[kIndex]['partial_genotypes'])
                positive_genos[base] = set(positive_genos[base])
                positive_genos[base] = set(positive_genos[base])
            optimal_starts = []
            selected_kmers = set()
            for i in range(0, num_bases):
                b1 = bases[i]
                if b1 not in self.kmer_scheme_data[pos]:
                    continue

                kpos_geno_counts = {}
                kpos_inf_scores = {}
                kpos_kmer_index = {}
                #loop through the kmer scheme dictionary to find the best kmer sequences
                #with the best starting positions
                for k in self.kmer_scheme_data[pos][b1]:
                    s = kmer_info[k]['aln_start']
                    if not s in kpos_geno_counts:
                        kpos_geno_counts[s] = {}
                        kpos_inf_scores[s] = 0
                        kpos_kmer_index[s] = []
                    kpos_kmer_index[s].append(k)
                    if len(kmer_rules[k]['positive_genotypes']) > 0:
                        kpos_inf_scores[s] += 1
                    counts = kmer_info[k]['genotype_counts']
                    for g in counts:
                        if not g in kpos_geno_counts[s]:
                            kpos_geno_counts[s][g] = 0
                        kpos_geno_counts[s][g] += counts[g]
                kpos_inf_scores = {k: v for k, v in
                                   sorted(kpos_inf_scores.items(), key=lambda item: item[1], reverse=True)}
                best_start = list(kpos_inf_scores.keys())[0]
                optimal_starts.append(best_start)
                selected_kmers = selected_kmers | set(kpos_kmer_index[best_start])
            #remove all but the best kmers found above
            for i in range(0, num_bases):
                b1 = bases[i]
                if b1 not in self.kmer_scheme_data[pos]:
                    continue
                self.kmer_scheme_data[pos][b1] = list(set(self.kmer_scheme_data[pos][b1]) & selected_kmers)


    def get_kseq_by_index(self,index):
        """
        A helper method to return the kmer sequence for the kmer index passed
        from the extracted kmers dictionary.

        Parameters
        ----------
        index - int
            The index for the desired kmer sequence

        Returns
        -------
        str
            The kmer sequence if it is found in the extracted kmer dictionary, but returns an empty string if not found
        """
        i = 0
        for kmer in self.extracted_kmers:
            if i == index:
                return kmer
            i+=1
        return ''

    def remove_scheme_pos(self,pos_to_remove):
        """
        A helper method to remove an entry in the kmer scheme data dictionary
        based on the position passed
        :param pos_to_remove: int - the index of the position to be removed
        """
        for pos in pos_to_remove:
            if pos in self.kmer_scheme_data:
                del(self.kmer_scheme_data[pos])


    def calc_consensus_seq(self):
        """
        A method to determine the consensus sequence for the fasta passed
        by looping through the reference sequence.

        Returns
        -------
        str
            The consensus sequence for the bases in the sequence fasta passed
        """
        consensus = []
        valid_bases = ['A','T','C','G']
        for pos in self.msa_base_counts:
            bases = []
            for b in valid_bases:
                count = self.msa_base_counts[pos][b]
                if count > 0:
                    bases.append(b)
            if len(bases) ==0:
                c = '-'
            else:
                bases = "".join(sorted(bases))
                c = IUPAC_LOOK_UP[bases]
            consensus.append(c)
        return "".join(consensus)

    def create_biohansel_kmers(self):
        """
        A method to create the kmer data dictionary for use with the BioHansel scheme.
        This method creates the positive and negative kmers for use in the BioHansel scheme,
        while adhering to the same kmer rules as the :meth:`cladeomatic.create.create_scheme` method.

        Returns
        -------
        dict
            A dictionary for the kmers to be used by the BioHansel scheme

        Notes
        -----
        Refer to https://github.com/phac-nml/biohansel for more thorough BioHansel documentation

        """
        kmers = {}
        valid_bases = ['A', 'T', 'C', 'G']
        consensus = list(self.calc_consensus_seq())
        ref = self.ref_sequence
        opt_kstart_pos = self.opt_kmer_start_positions
        msa_len = len(consensus)
        klen = self.kmer_len
        for pos in self.target_positions:
            #skip any missing kmers
            if pos in self.positions_missing_kmer:
                continue
            kmers[pos] = {
                'A':{'positive':'','negative':'','start':-1,'end':-1,'genotype':''},
                'T':{'positive':'','negative':'','start':-1,'end':-1,'genotype':''},
                'C':{'positive':'','negative':'','start':-1,'end':-1,'genotype':''},
                'G':{'positive':'','negative':'','start':-1,'end':-1,'genotype':''},
            }
            ref_base = ref[pos]
            cons_base = consensus[pos]
            bases = []
            for b in valid_bases:
                count = self.msa_base_counts[pos][b]
                if count > 0:
                    bases.append(b)
            alt_bases = set(bases) - set(ref_base)
            num_alt_bases = len(alt_bases)
            kStart = pos - klen + 1
            if kStart < 0:
                kStart = 0
            kEnd = pos + klen + 1
            if kEnd >= msa_len:
                kEnd = msa_len - 1
            pos_ovl_variant = set(range(kStart,kEnd) ) & set(opt_kstart_pos)
            if len(pos_ovl_variant) == 0:
                pos_ovl_variant = set(kStart)
            pos_ovl_variant = list(pos_ovl_variant)
            if len(pos_ovl_variant) < num_alt_bases:
                for i in range(0,num_alt_bases-len(pos_ovl_variant)):
                    kStart = pos_ovl_variant[-1] + 1
                    pos_ovl_variant.append(kStart)

            if num_alt_bases == 0:
                kStart = pos_ovl_variant[0]
                kEnd = kStart + klen + 1
                ref_kmer = consensus[kStart:kEnd]
                kmers[pos][base]['positive'] = ''.join(ref_kmer)
                kmers[pos][base]['start'] = kStart
                kmers[pos][base]['end'] = kEnd
                continue
            alt_bases = list(alt_bases)
            for i in range(0,len(alt_bases)):
                kStart = pos_ovl_variant[i]
                kEnd = kStart + klen + 1
                base = alt_bases[i]
                consensus[pos] = base
                alt_kmer = ''.join(consensus[kStart:kEnd])
                out_bases = set(bases) - set(base)
                iupac_key = []
                for b in valid_bases:
                    if b in out_bases:
                        iupac_key.append(b)
                iupac_key = ''.join(iupac_key)
                consensus[pos] = IUPAC_LOOK_UP[iupac_key]
                ref_kmer = consensus[kStart:kEnd]
                kmers[pos][base]['positive'] = alt_kmer
                kmers[pos][base]['negative'] = ''.join(ref_kmer)
                kmers[pos][base]['start'] = kStart
                kmers[pos][base]['end'] = kEnd

            consensus[pos] = cons_base

        return kmers





