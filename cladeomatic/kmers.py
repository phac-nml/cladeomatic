import time

from Bio import SeqIO

from cladeomatic.utils.kmerSearch import SeqSearchController


class kmer_worker:
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


    def __init__(self, ref_sequence, msa_file, result_dir, prefix, klen, genotype_map, max_ambig=0, min_perc=1,
                 target_positions=[], num_threads=1):
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
        self.workflow()
        return

    def workflow(self):
        print(self.target_positions)
        stime = time.time()
        self.init_msa_base_counts()
        print(stime - time.time())
        stime = time.time()
        self.pop_msa_base_counts()
        print(stime - time.time())
        stime = time.time()
        self.find_variant_positions()
        print(stime - time.time())
        if len(self.target_positions) == 0:
            self.target_positions = self.variant_positions
        stime = time.time()
        self.get_genotype_snp_states()
        print(stime - time.time())
        stime = time.time()
        self.get_optimal_kmer_position()
        print(stime - time.time())
        stime = time.time()
        self.extract_kmers()
        print(stime - time.time())
        stime = time.time()
        self.populate_int_base_kmer_lookup()
        print(stime - time.time())
        stime = time.time()
        self.perform_kmer_search()
        print(stime - time.time())
        stime = time.time()
        self.process_kmer_results()
        print(stime - time.time())
        stime = time.time()
        self.init_kmer_scheme_data()
        print(stime - time.time())
        stime = time.time()
        self.populate_kmer_scheme_data()
        print(stime - time.time())
        stime = time.time()

        self.confirm_kmer_specificity()
        print(stime - time.time())

        # Reinitialize with invalid kmers removed
        stime = time.time()
        self.init_kmer_scheme_data()
        print(stime - time.time())
        stime = time.time()
        self.populate_kmer_scheme_data()
        print(stime - time.time())
        stime = time.time()
        self.remove_empty_base_states()
        print(stime - time.time())
        stime = time.time()
        self.construct_ruleset()
        print(stime - time.time())
        stime = time.time()
        self.remove_redundant_kmers()
        print(stime - time.time())
        stime = time.time()
        self.find_invalid_kmers()
        print(stime - time.time())
        stime = time.time()
        self.remove_invalid_kmers_from_scheme()
        print(stime - time.time())
        stime = time.time()
        self.positions_missing_kmer = self.get_pos_without_kmer()
        print(stime - time.time())

    def init_msa_base_counts(self):
        for i in range(0, self.ref_len):
            self.msa_base_counts[i] = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, '-': 0}

    def pop_msa_base_counts(self):
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
        for pos in self.msa_base_counts:
            count_bases = 0
            for base in self.valid_bases:
                if self.msa_base_counts[pos][base] > 0:
                    count_bases += 1
            if count_bases > 1:
                self.variant_positions.append(pos)

    def get_optimal_kmer_position(self):
        msa_len = self.ref_len
        klen = self.kmer_len
        local_kmer_start_positions = set()
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
            for i in interval:
                variable_sites = 0
                num_missing = 0
                for k in range(i, i + klen + 1):
                    if k >= end:
                        break
                    num_bases = 0
                    if self.msa_base_counts[k]['-'] > 0:
                        num_missing += 1

                    for base in self.valid_bases:
                        if self.msa_base_counts[k][base] > 0:
                            num_bases += 1
                    if num_bases > 1:
                        variable_sites += 1

                if variable_sites <= best_count_variable_sites:
                    if num_missing < best_count_missing_sites:
                        best_index = i
                        best_count_variable_sites = variable_sites
                        best_count_missing_sites = num_missing
                local_kmer_start_positions.add(best_index)
        self.opt_kmer_start_positions = sorted(list(local_kmer_start_positions))

    def extract_kmers(self):
        msa_len = self.ref_len
        klen = self.kmer_len
        selected_kmers = {}
        canonical_pos = set(self.target_positions)
        with open(self.msa_fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper()
                for pos in self.variant_positions:
                    base = seq[pos]
                    if base not in self.valid_bases:
                        continue
                    s = pos - klen + 1
                    if s < 0:
                        s = 0
                    start_range = (s, pos)
                    for start in start_range:
                        end = start + klen
                        kseq = seq[start:end + 1]
                        min_ambig = kseq.count('N')
                        no_gap_len = len(kseq.replace("-", ""))
                        while no_gap_len > klen and end > pos and kseq.count('N') <= min_ambig:
                            end -= 1
                            kseq = seq[start:end + 1]
                            no_gap_len = len(kseq.replace("-", ""))
                            if kseq.count('N') <= min_ambig:
                                min_ambig = kseq.count('N')

                        while no_gap_len < klen and end < msa_len and kseq.count('N') <= min_ambig:
                            end += 1
                            kseq = seq[start:end + 1]
                            no_gap_len = len(kseq.replace("-", ""))
                            if kseq.count('N') <= min_ambig:
                                min_ambig = kseq.count('N')

                        while no_gap_len < klen and start > 0 and kseq.count('N') <= min_ambig:
                            start -= 1
                            kseq = seq[start:end + 1]
                            no_gap_len = len(kseq.replace("-", ""))
                            if kseq.count('N') <= min_ambig:
                                min_ambig = kseq.count('N')

                        kseq = seq[start:end + 1].replace('-', '')
                        if kseq.count('N') > self.max_ambig or len(kseq) != klen:
                            continue
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
        seqkmers = {}
        index = 0
        for kmer in self.extracted_kmers:
            seqkmers[index] = kmer
            index += 1
        self.kmer_search_files = SeqSearchController(seqkmers, self.msa_fasta_file, self.result_dir, self.prefix,
                                                     self.num_threads)

    def process_kmer_results(self):
        invalid_kmers = set()
        genotype_mapping = self.genotype_membership
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
                if count > 1:
                    invalid_kmers.add(kIndex)
                if not genoytpe in self.int_base_kmer_lookup[kIndex]['genotype_counts']:
                    self.int_base_kmer_lookup[kIndex]['genotype_counts'][genoytpe] = 0
                self.int_base_kmer_lookup[kIndex]['genotype_counts'][genoytpe] += 1
            self.flag_invalid_kmers(invalid_kmers)

    def flag_invalid_kmers(self, invalid_kmers):
        kIndex = 0
        for kmer in self.extracted_kmers:
            if kIndex in invalid_kmers:
                self.extracted_kmers[kmer]['is_valid'] = False
            kIndex += 1

    def init_kmer_scheme_data(self):
        kmer_data = {}
        for pos in self.target_positions:
            kmer_data[pos] = {
                'A': [], 'T': [], 'C': [], 'G': [],
            }
        self.kmer_scheme_data = kmer_data

    def populate_kmer_scheme_data(self):
        scheme_data = self.kmer_scheme_data
        kmer_data = self.extracted_kmers
        index = 0
        invalid_kmers = set()
        for kmer in kmer_data:
            is_valid = kmer_data[kmer]['is_valid']
            ovl_pos = kmer_data[kmer]['target_positions']
            if not is_valid:
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
        genotype_map = self.genotype_membership
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

    def get_genotype_snp_states(self):
        with open(self.msa_fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                id = str(record.id).split("~")[0]
                if id not in self.genotype_membership:
                    continue
                seq = str(record.seq).upper()
                genotype = self.genotype_membership[id]
                if not genotype in self.genotype_snp_associations:
                    self.genotype_snp_associations[genotype] = {}
                for pos in self.target_positions:
                    if not pos in self.genotype_snp_associations[genotype]:
                        self.genotype_snp_associations[genotype][pos] = set()
                    base = seq[pos]
                    if not base in self.valid_bases:
                        continue
                    self.genotype_snp_associations[genotype][pos].add(base)
            handle.close()

    def confirm_kmer_specificity(self):
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
        index = 0
        for kmer in self.extracted_kmers:
            if self.extracted_kmers[kmer]['is_valid'] == False:
                self.invalid_kmer_indexes.add(index)
            index += 1

    def populate_int_base_kmer_lookup(self):
        index = 0
        for kmer in self.extracted_kmers:
            self.int_base_kmer_lookup[index] = self.extracted_kmers[kmer]
            index += 1

    def find_ovl_kmers(self):
        variant_postions = self.target_positions
        selected_kmers = self.kmer_scheme_data
        num_var_pos = len(variant_postions)
        shared_kmers = set()
        pot_pos_ovl = {}
        for i in range(0, len(variant_postions)):
            var_pos_1 = variant_postions[i]
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
        for index in self.kmer_scheme_data:
            for base in self.kmer_scheme_data[index]:
                ovl = set(self.kmer_scheme_data[index][base]) - self.invalid_kmer_indexes
                self.kmer_scheme_data[index][base] = ovl

    def remove_empty_base_states(self):
        for index in self.kmer_scheme_data:
            for base in self.valid_bases:
                if base not in self.kmer_scheme_data[index]:
                    continue
                if len(self.kmer_scheme_data[index][base]) == 0:
                    del (self.kmer_scheme_data[index][base])

    def remove_redundant_kmers(self):
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

            for i in range(0, num_bases):
                b1 = bases[i]
                if b1 not in self.kmer_scheme_data[pos]:
                    continue

                kpos_geno_counts = {}
                kpos_inf_scores = {}
                kpos_kmer_index = {}

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
                self.kmer_scheme_data[pos][b1] = list(kpos_kmer_index[best_start])

    def get_kseq_by_index(self,index):
        i = 0
        for kmer in self.extracted_kmers:
            if i == index:
                return kmer
            i+=1
        return ''

    def remove_scheme_pos(self,pos_to_remove):
        for pos in pos_to_remove:
            if pos in self.kmer_scheme_data:
                del(self.kmer_scheme_data[pos])







