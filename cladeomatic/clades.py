import pandas as pd
import numpy as np
from cladeomatic.snps import snp_search_controller
from cladeomatic.utils import calc_ARI, calc_AMI, calc_shanon_entropy, fisher_exact
from scipy.signal import find_peaks
from scipy.spatial import distance
from scipy.cluster import hierarchy
from collections import OrderedDict
import math

class clade_worker:
    vcf_file = None
    metadata_dict = {}
    distance_matrix_file = None
    group_data = {}
    delim = '.'
    perform_compression = True
    num_threads = 1
    min_member_count = 1
    min_snp_count = 1
    max_states = 6
    min_inter_clade_dist = 1
    ref_seq = ''
    method = ''

    #Derived
    snp_data = {}
    clade_data = {}
    raw_genotypes = {}
    supported_genotypes = {}
    selected_genotypes = {}
    node_list = []
    node_counts = {}
    bifurcating_nodes = {}
    terminal_nodes = set()
    node_bins = []
    sample_distance_histo = {}
    sample_linkage_clusters = {}
    dist_thresholds = []
    dist_ranges = []
    variant_positions = []
    is_root_valid = True
    selected_positions = []


    def __init__(self,vcf,metadata_dict,dist_mat_file,groups,ref_seq,perform_compression=True,delim='.',min_snp_count=1,max_states=6,min_members=1,min_inter_clade_dist=1,num_threads=1,method='average'):
        self.vcf_file = vcf
        self.metadata_dict = metadata_dict
        self.ref_seq = ref_seq
        self.distance_matrix_file = dist_mat_file
        self.group_data = groups
        self.num_threads = num_threads
        self.min_member_count = min_members
        self.delim = delim
        self.max_states = max_states
        self.min_snp_count = min_snp_count
        self.min_inter_clade_dist = min_inter_clade_dist
        self.perform_compression = perform_compression
        self.method = method
        self.workflow()
        return

    def workflow(self):
        self.raw_genotypes = self.generate_genotypes()
        self.snp_data = snp_search_controller(self.group_data, self.vcf_file, self.num_threads)
        self.variant_positions = self.get_variant_positions()
        self.snp_based_filter()
        self.get_all_nodes()
        self.get_node_member_counts()
        self.init_clade_data()
        self.populate_clade_data()
        self.calc_node_distances()
        self.check_nodes()
        valid_nodes = self.get_valid_nodes()
        self.set_valid_nodes(valid_nodes)
        self.supported_genotypes = self.generate_genotypes()
        self.calc_node_associations_groups()


        if self.perform_compression:
            self.set_invalid_nodes(self.get_close_nodes())
            self.get_bifurcating_nodes()
            self.distance_node_binning()
            valid_nodes = self.get_valid_nodes()
            self.set_valid_nodes(valid_nodes)
            self.terminal_nodes = self.get_terminal_nodes()
            self.set_invalid_nodes( self.get_close_nodes())
            self.compress_heirarchy()
            self.fix_root()
            valid_nodes = self.get_valid_nodes()
            self.set_valid_nodes(valid_nodes)
            self.dist_based_nomenclature()

        self.selected_genotypes = self.generate_genotypes()
        self.selected_positions = self.get_selected_positions()


    def get_selected_positions(self):
        selected_positions = set()
        for node_id in self.clade_data:
            if self.clade_data[node_id]['is_valid']:
                selected_positions = selected_positions | set(self.clade_data[node_id]['pos'])
        return sorted(list(selected_positions))

    def set_valid_nodes(self,valid_nodes):
        self.group_data['valid_nodes'] = set(valid_nodes)

    def set_invalid_nodes(self,invalid_nodes):
        for node_id in invalid_nodes:
            self.clade_data[node_id]['is_valid'] = False
            self.clade_data[node_id]['is_selected'] = False

    def get_all_nodes(self):
        self.node_list = sorted(list(self.group_data['valid_nodes']))

    def get_valid_nodes(self):
        valid_nodes = set()
        for node_id in self.clade_data:
            if self.clade_data[node_id]['is_valid']:
                valid_nodes.add(node_id)
        return valid_nodes

    def get_terminal_nodes(self):
        genotypes = self.generate_genotypes()
        terminal_nodes = set()
        for sample_id in genotypes:
            terminal_nodes.add(genotypes[sample_id].split(self.delim)[-1])

        return terminal_nodes

    def get_selected_nodes(self):
        valid_nodes = set()
        for node_id in self.clade_data:
            if self.clade_data[node_id]['is_selected']:
                valid_nodes.add(node_id)
        return valid_nodes

    def get_node_member_counts(self):
        counts = {}
        for sample_id in self.group_data:
            for n in self.group_data[sample_id]:
                if n not in counts:
                    counts[n] = 0
                counts[n]+=1
        self.node_counts = counts

    def init_clade_data(self):
        for node_id in self.node_list:
            self.clade_data[node_id] = {
                'pos': [],
                'bases': [],
                'min_dist': 0,
                'max_dist': 0,
                'ave_dist': 0,
                'entropies': {},
                'fisher': {},
                'ari': {},
                'ami': {},
                'num_members':self.node_counts[node_id],
                'is_valid':True,
                'is_selected':True,
                'total_within_clade_dist':0,
                'num_comparisons_within_clade_dist': 0,
                'ave_within_clade_dist': 0,
                'closest_clade_id':None,
                'closest_clade_dist':0,
                'closest_sample_id': None,
                'closest_sample_dist': 0,
                'clade_sample_id':None
            }

    def populate_clade_data(self):
        for chrom in self.snp_data:
            for pos in self.snp_data[chrom]:
                for base in self.snp_data[chrom][pos]:
                    if self.snp_data[chrom][pos][base]['is_canonical']:
                        clade_id = self.snp_data[chrom][pos][base]['clade_id'].split('.')[-1]
                        self.clade_data[clade_id]['pos'].append(pos)
                        self.clade_data[clade_id]['bases'].append(base)

    def check_nodes(self):
        for node_id in self.clade_data:
            if self.clade_data[node_id]['num_members'] < self.min_member_count:
                self.clade_data[node_id]['is_valid'] = False
                self.clade_data[node_id]['is_selected'] = False
            if len(self.clade_data[node_id]['pos']) < self.min_snp_count:
                self.clade_data[node_id]['is_valid'] = False
                self.clade_data[node_id]['is_selected'] = False

    def generate_genotypes(self):
        '''

        Parameters
        ----------
        group_data
        delim

        Returns
        -------

        '''
        if self.delim == None:
            self.delim = '.'

        genotypes = {}
        for id in self.group_data['sample_map']:
            sample_id = self.group_data['sample_map'][id]['sample_id']
            genotype = [str(x) for x in self.group_data['sample_map'][id]['genotype'].split(self.delim)]
            filt_genotype = []
            for i in range(0, len(genotype)):
                node_id = genotype[i]
                if node_id in self.group_data['valid_nodes']:
                    filt_genotype.append(node_id)

            genotypes[sample_id] = "{}".format(self.delim).join(filt_genotype)
        return genotypes

    def snp_based_filter(self):
        '''
        Accepts a snp_data dict and filters the dictionary to return only those snps which are present in >=min_count samples
        and have no more than max_state bases at a given position
        :param snp_data: dict() {[chrom_id] : {[position]: {[base]: :dict()}}}
        :param min_count: int
        :param max_states: int
        :return: dict
        '''
        filtered = {}
        for chrom in self.snp_data:
            filtered[chrom] = {}
            for pos in self.snp_data[chrom]:
                for base in self.snp_data[chrom][pos]:
                    if self.snp_data[chrom][pos][base]['is_canonical']:
                        if self.snp_data[chrom][pos][base]['num_members'] >= self.min_member_count:
                            if not pos in filtered[chrom]:
                                filtered[chrom][pos] = {}
                            filtered[chrom][pos][base] = self.snp_data[chrom][pos][base]
                if pos in filtered[chrom] and len(filtered[chrom][pos]) > self.max_states:
                    del (filtered[chrom][pos])
        self.snp_data = filtered

    def get_bifurcating_nodes(self):
        nodes = {}
        valid_nodes = self.get_valid_nodes()
        genotypes = self.generate_genotypes()
        for sample_id in genotypes:
            genotype = genotypes[sample_id].split(self.delim)
            num_parent_nodes = len(genotype) - 1
            for i in range(0, num_parent_nodes):
                if i == num_parent_nodes:
                    continue
                parent_node_id = genotype[i]
                child_node_id = genotype[i+1]
                if not parent_node_id in nodes:
                    nodes[parent_node_id] = set()
                nodes[parent_node_id].add(child_node_id)

        for node_id in nodes:
            if len(nodes[node_id]) > 1:
                self.bifurcating_nodes[node_id] = nodes[node_id]

    def genotype_lookup(self,sample_list):
        lookup = {}
        genotypes = self.generate_genotypes()
        d = self.delim
        for sample_id in sample_list:
            lookup[sample_id] = genotypes[sample_id].split(d)

        return lookup

    def calc_node_distances(self):
        fh = open(self.distance_matrix_file,'r')
        header = next(fh).rstrip().split("\t")
        num_columns = len(header)
        sample_list = header[1:]
        sample_lookup = self.genotype_lookup(sample_list)
        histo = {}
        index = 1
        for line in fh:
            line = line.rstrip().split("\t")
            sample_id_1 = line[0]
            nodes_1 = set(sample_lookup[sample_id_1])
            for i in range(index,num_columns):
                sample_id_2 = header[i]
                if sample_id_1 == sample_id_2:
                    continue
                nodes_2 = set(sample_lookup[sample_id_2])
                value = int(line[i])
                if not value in histo:
                    histo[value] = 0;
                histo[value] +=1
                ovl = nodes_1 & nodes_2
                for node_id in ovl:
                    self.clade_data[node_id]['total_within_clade_dist'] += value
                    self.clade_data[node_id]['num_comparisons_within_clade_dist'] += 1
                no_ovl_1 =  nodes_1 - nodes_2
                for node_id in no_ovl_1:
                    if self.clade_data[node_id]['closest_sample_id'] is None or self.clade_data[node_id]['closest_sample_dist'] < value:
                        self.clade_data[node_id]['clade_sample_id'] = sample_id_1
                        self.clade_data[node_id]['closest_sample_id'] = sample_id_2
                        self.clade_data[node_id]['closest_sample_dist'] = value

                no_ovl_2 =  nodes_2 - nodes_1
                for node_id in no_ovl_2:
                    if self.clade_data[node_id]['closest_sample_id'] is None or self.clade_data[node_id][
                        'closest_sample_dist'] < value:
                        self.clade_data[node_id]['clade_sample_id'] = sample_id_2
                        self.clade_data[node_id]['closest_sample_id'] = sample_id_1
                        self.clade_data[node_id]['closest_sample_dist'] = value

            index+=1

        fh.close()
        for node_id in self.clade_data:
            if self.clade_data[node_id]['num_comparisons_within_clade_dist'] > 0:
                ave = self.clade_data[node_id]['total_within_clade_dist'] / self.clade_data[node_id]['num_comparisons_within_clade_dist']
                self.clade_data[node_id]['ave_within_clade_dist'] = ave
            sample_id_1 = self.clade_data[node_id]['clade_sample_id']
            sample_id_2 = self.clade_data[node_id]['closest_sample_id']
            if sample_id_1 is None:
                continue
            nodes_1 = sample_lookup[sample_id_1]
            nodes_2 = sample_lookup[sample_id_2]
            num_nodes = len(nodes_1)
            if len(nodes_2) < num_nodes:
                num_nodes = len(nodes_2)
            shared_parent = None
            for i in range(0,num_nodes):
                n1 = nodes_1[i]
                n2 = nodes_2[i]
                if n1 == n2:
                    shared_parent = n1
            self.clade_data[node_id]['closest_clade_id'] = shared_parent
            self.clade_data[node_id]['closest_clade_dist'] = self.clade_data[node_id]['closest_sample_dist']



        self.distance_histo = histo

    def get_clade_distances(self):
        distances = []
        for node_id in self.clade_data:
            distances.append(self.clade_data[node_id]['ave_within_clade_dist'])
        return distances

    def get_inter_clade_distances(self):
        distances = []
        for node_id in self.clade_data:
            distances.append(self.clade_data[node_id]['closest_clade_dist'])
        return distances

    def get_close_nodes(self):
        invalid_nodes = set()
        for node_id in self.clade_data:
            dist = self.clade_data[node_id]['closest_clade_dist']
            if dist < self.min_inter_clade_dist:
                invalid_nodes.add(node_id)
                if self.clade_data[node_id]['closest_clade_id'] is not None:
                    invalid_nodes.add(self.clade_data[node_id]['closest_clade_id'])
        return invalid_nodes

    def partition_distances(self):
        num_partitions = 10
        distances = sorted(self.get_clade_distances())
        num_distances = len(distances)
        target_number_samples = int(num_distances / num_partitions)

        bins = {}
        for i in range(0,num_partitions):
            bins[i] = []
        bin_index = 0

        for i in range(0,num_distances):
            d = distances[i]

            if len(bins[bin_index]) > target_number_samples:
                if d > distances[i-1]:
                    bin_index+=1
            bins[bin_index].append(d)

        partitions = []
        for i in bins:
            if len(bins[i]) > 0:
                partitions.append(bins[i][-1])

        return partitions

    def distance_node_binning(self):
        partition_dists = self.partition_distances()
        num_bins = len(partition_dists)
        bins = {}
        for i in range(0,num_bins):
            bins[i] = set()
        bin_indexes = range(0,num_bins )
        for node_id in self.node_list:
            d = self.clade_data[node_id]['ave_within_clade_dist']
            for i in bin_indexes:
                t = partition_dists[i]
                if d <= t:
                    bins[i].add(node_id)
        self.node_bins = bins

    def compress_heirarchy(self):
        node_bins = self.node_bins
        genotypes = self.generate_genotypes()
        invalid_nodes = set()
        for sample_id in genotypes:
            genotype = genotypes[sample_id].split(self.delim)
            num_nodes = len(genotype)
            if num_nodes < 2:
                continue
            for i in range(0,num_nodes-1):
                n1 = genotype[i]
                n2 = genotype[i+1]
                if n2 in self.bifurcating_nodes:
                    continue
                for bin_id in node_bins:
                    if n1 in node_bins[bin_id] and n2 in node_bins[bin_id]:
                        invalid_nodes.add(n2)
                        break

        self.set_invalid_nodes(invalid_nodes)
        self.set_valid_nodes(self.get_valid_nodes() - invalid_nodes)

    def get_variant_positions(self):
        variant_positions = []
        for chrom in self.snp_data:
            for pos in self.snp_data[chrom]:
                variant_positions.append(pos)

        return variant_positions

    def get_conserved_ranges(self):
        vpos = self.variant_positions
        num_pos = len(vpos)
        if vpos[0] > 0:
            ranges = [[0,vpos[0]-1]]
        ranges = []
        for i in range(0,num_pos-1):
            p1 = vpos[i]
            for k in range(i+1,num_pos):
                p2 = vpos[k]
                if p2 - p1 == 0:
                    break
                ranges.append([p1+1,p2-1])
        ref_len = len(self.ref_seq)
        if vpos[-1] < ref_len -1:
            ranges.append([vpos[-1] + 1, ref_len - 1])

        return ranges

    def get_largest_conserved_pos(self):
        ranges = self.get_conserved_ranges()
        if len(ranges) == 0:
            self.is_root_valid = False
            return
        gap_len = 0
        index = 0
        for i in range(0,len(ranges)):
            start = ranges[i][0]
            end = ranges[i][1]
            l = end - start
            if l > gap_len:
                index = i
                gap_len = l

        if gap_len == 1:
            pos = ranges[index][0]
        else:
            pos = ranges[index][0] + int(gap_len/2)
        return pos

    def fix_root(self):
        pos = self.get_largest_conserved_pos()
        base = self.ref_seq[pos]
        self.clade_data['0']['pos'].append(pos)
        self.clade_data['0']['bases'].append(base)
        self.clade_data['0']['is_valid'] = True
        self.clade_data['0']['is_selected'] = True

    def perform_clustering(self):
        data = pd.read_csv(self.distance_matrix_file, sep='\t', header=0, index_col=0)
        data_labels = data.columns.values
        distance_matrix = distance.squareform(data.values)
        Z = hierarchy.linkage(distance_matrix, method=self.method)
        cluster_membership = {}
        for id in data_labels:
            cluster_membership[id] = []
        clust_number = 1
        for dist in self.dist_thresholds:
            clusters = hierarchy.fcluster(Z, dist, criterion='distance')
            ids = set(clusters)
            id_map = {}
            for id in ids:
                id_map[id] = clust_number
                clust_number+=1
            index = 0
            for id in data_labels:
                cluster_membership[id].append(id_map[clusters[index]])
                index+=1
        self.sample_linkage_clusters = cluster_membership

    def create_cluster_membership_lookup(self):
        nodes = {}
        sample_map = self.group_data['sample_map']
        for index in sample_map:
            sample_id = sample_map[index]['sample_id']
            for node_id in self.sample_linkage_clusters[sample_id]:
                if not node_id in nodes:
                    nodes[node_id] = set()
                nodes[node_id].add(index)
        return nodes

    def fit_clusters_to_groupings(self):
        cluster_to_node_mapping = {}
        sl_cluster_membership = self.create_cluster_membership_lookup()
        for sl_clust_id in sl_cluster_membership:
            clust_members = sl_cluster_membership[sl_clust_id]
            num_members = len(clust_members)
            best_node_id = '0'
            best_node_nonmembers = len(self.group_data['membership'])
            for node_id in self.group_data['membership']:
                node_members = self.group_data['membership'][node_id]
                missing = clust_members - node_members
                if len(missing) > 0:
                    continue
                nonmembers = node_members - clust_members
                count_nonmembers = len(nonmembers)
                if count_nonmembers < best_node_nonmembers:
                    best_node_id = node_id
                    best_node_nonmembers = count_nonmembers
                if count_nonmembers == 0:
                    break
            cluster_to_node_mapping[sl_clust_id] = best_node_id.split(self.delim)[-1]
        return cluster_to_node_mapping

    def identify_dist_ranges(self):
        histo = OrderedDict(sorted(self.distance_histo.items()))
        histo_keys = list(histo.keys())
        max_value = histo_keys[-1]
        x = list(range(0,max_value+2))
        y = []
        for i in x:
            value = 0
            if i in histo:
                value = histo[i]
            y.append(value)
        troughs, properties = self.find_dist_troughs(y)
        self.dist_thresholds = sorted(list(troughs),reverse=True)
        troughs = sorted(list(set(troughs)))
        c = self.as_range(troughs)
        ranges = []
        for idx, tuple in enumerate(c):
            ranges.append((troughs[tuple[0]], troughs[tuple[1]]))

        return ranges

    def dist_based_nomenclature(self):
        self.dist_ranges = self.identify_dist_ranges()
        self.perform_clustering()
        clust_to_node_mapping = self.fit_clusters_to_groupings()
        sample_genotypes = {}
        for sample_id in self.sample_linkage_clusters:
            clust_geno = self.sample_linkage_clusters[sample_id]
            node_geno = []
            for clust_id in  clust_geno:
                node_id = clust_to_node_mapping[clust_id]
                node_geno.append(node_id)
            if len(node_geno) > 1:
                filt = [node_geno[0]]
                for i in range(1,len(node_geno)-1):
                    if node_geno[i] != node_geno[i-1]:
                        filt.append(node_geno[i])
            else:
                filt = node_geno
            sample_genotypes[sample_id] = self.delim.join(filt)

    def find_dist_troughs(self,values):
        '''
        :param values: list
        :return:
        '''
        return (find_peaks(-np.array(values)))

    def as_range(self,values):
        '''
        Parameters
        ----------
        values list of integers
        Returns list of tuples with the start and end index of each range
        -------
        '''
        values = sorted(values)
        num_val = len(values)
        if num_val == 0:
            return []
        if num_val == 1:
            return [(0, 0)]
        ranges = []
        i = 0
        while i < num_val - 1:

            s = i
            for k in range(i + 1, num_val):
                if values[k] - values[i] != 1:
                    i += 1
                    break
                i += 1
            if values[k] - values[i] > 1:
                ranges.append((s, k - 1))
            else:
                if k == num_val - 1:
                    ranges.append((s, k))
                else:
                    ranges.append((s, k - 1))
        return ranges

    def remove_snps(self,positions):
        positions = set(positions)
        for node_id in self.clade_data:
            ovl = set(self.clade_data[node_id]['pos']) & positions
            if len(ovl) == 0:
                continue
            num_pos = len(self.clade_data[node_id]['pos'])
            p = []
            b = []
            for i in range(0,num_pos):
                pos = self.clade_data[node_id]['pos'][i]
                if pos in ovl:
                    continue
                p.append(pos)
                b.append(self.clade_data[node_id]['bases'][i])
            self.clade_data[node_id]['pos'] = p
            self.clade_data[node_id]['bases'] = b
            if len(p) == 0:
                self.clade_data[node_id]['is_valid'] = False
                self.clade_data[node_id]['is_selected'] = False
        self.set_valid_nodes(self.get_valid_nodes())
        self.selected_genotypes = self.generate_genotypes()

    def calc_node_associations_groups(self):
        '''
        :param metadata:
        :param clade_data:
        :param ete_tree_obj:
        :return:
        '''

        samples = set(self.metadata_dict.keys())
        num_samples = len(samples)
        for clade_id in self.group_data['membership']:
            if not clade_id in self.clade_data:
                continue
            in_members = self.group_data['membership'][clade_id]
            features = {}
            genotype_assignments = []
            metadata_labels = {}
            ftest = {}
            metadata_counts = {}
            for sample_id in samples:
                if sample_id in in_members:
                    genotype_assignments.append(1)
                else:
                    genotype_assignments.append(0)
                for field_id in self.metadata_dict[sample_id]:
                    value = self.metadata_dict[sample_id][field_id]
                    if not field_id in metadata_labels:
                        metadata_counts[field_id] = {}
                        metadata_labels[field_id] = []
                        ftest[field_id] = {}
                    if not value in ftest[field_id]:
                        metadata_counts[field_id][value] = 0
                        ftest[field_id][value] = {
                            'pos-pos': set(),
                            'pos-neg': set(),
                            'neg-pos': set(),
                            'neg-neg': set()
                        }
                    metadata_counts[field_id][value] += 1
                    if sample_id in in_members:
                        ftest[field_id][value]['pos-pos'].add(sample_id)
                    else:
                        ftest[field_id][value]['neg-pos'].add(sample_id)

                    metadata_labels[field_id].append(value)
                    if not field_id in features:
                        features[field_id] = {}

                    if not value in features[field_id]:
                        features[field_id][value] = 0
                    if sample_id in in_members:
                        features[field_id][value] += 1

            for field_id in metadata_labels:
                category_1 = []
                category_2 = []
                for idx, value in enumerate(metadata_labels[field_id]):
                    if isinstance(value, float):
                        if math.isnan(value):
                            continue
                    value = str(value)
                    if len(value) == 0 or value == 'nan':
                        continue
                    category_1.append(genotype_assignments[idx])
                    category_2.append(metadata_counts[field_id][value])

                self.clade_data[clade_id]['ari'][field_id] = calc_ARI(category_1, category_2)
                self.clade_data[clade_id]['ami'][field_id] = calc_AMI(category_1, category_2)
                self.clade_data[clade_id]['entropies'][field_id] = calc_shanon_entropy(category_2)
                self.clade_data[clade_id]['fisher'][field_id] = {}

                for value in ftest[field_id]:
                    ftest[field_id][value]['neg-neg'] = (
                            ftest[field_id][value]['pos-pos'] | ftest[field_id][value]['neg-pos'])
                    ftest[field_id][value]['pos-neg'] = in_members - ftest[field_id][value]['neg-neg']
                    table = [
                        [len(ftest[field_id][value]['pos-pos']),
                         len(ftest[field_id][value]['neg-neg'] | ftest[field_id][value]['pos-neg'])],
                        [len(ftest[field_id][value]['neg-pos'] | ftest[field_id][value]['pos-pos']),
                         len(ftest[field_id][value]['neg-neg'] | ftest[field_id][value]['neg-pos'])]
                    ]
                    oddsr, p = fisher_exact(table, alternative='greater')
                    self.clade_data[clade_id]['fisher'][field_id][value] = {'oddsr': oddsr, 'p': p}



























