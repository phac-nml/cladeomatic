import random
from collections import OrderedDict
from cladeomatic.utils.vcfhelper import vcfReader
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.signal import find_peaks
from scipy.spatial import distance
from scipy.stats import spearmanr, pearsonr

from cladeomatic.snps import snp_search_controller
from cladeomatic.utils import fisher_exact

from deprecated import deprecated

class clade_worker:
    """
    The clade worker class provides methods to take the user
    provided input and cluster the data, create and compress clades,
    and further summarize the SNP variants within the samples.
    """
    #initialization of variables
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
    rcor_thresh = 0
    max_snp_count = 1
    max_snp_resolution_thresh = 10
    min_perc = 1
    mode = None

    #Derived
    snp_data = {}
    selected_pos_ref_states = {}
    genotype_snp_data = {}
    num_positions = 0
    num_valid_positions = 0
    num_canonical_positions = 0
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
    sample_dist_based_nomenclature = {}
    cluster_to_node_mapping = {}
    genotype_snp_rules = {}


    def __init__(self,vcf,metadata_dict,dist_mat_file,groups,ref_seq,mode,perform_compression=True,delim='.',
                 min_snp_count=1,max_snps=-1,max_states=6,min_members=1,min_inter_clade_dist=1,num_threads=1,
                 max_snp_resolution_thresh=0,method='average',rcor_thresh=0.4,min_perc=1):
        """
        The instantiation of the :class:`clade_worker` class.

        Parameters
        ----------
        vcf : string
            The file path to the VCF file for clade processing
        metadata_dict : dict
            The dictionary of the metadata for the samples
        dist_mat_file : str
            The file path to the distance matrix file created by clade-o-matic
        groups : dict
            The dictionary containing the grouping of the data
        ref_seq : str
            The reference sequence
        mode : str
            A string to denote either 'tree' or 'group' mode for processing
        perform_compression : bool
            A boolean to denote if tree compression should occur.  Default is True.
        delim : str
            A string to denote the file/sample delimiter.  Default is '.'
        min_snp_count : int
            The minimum number of unique snps required for a clade to be valid.  Default is 1.
        max_snps : int
            The maximum number of snps required for the definition of a unique genotype.  Default is -1.
        max_states : int
            The maximum number of states for a position [A,T,C,G,N,-].  Default is 6.
        min_members : int
            The minimum number of members required for a clade to be considered valid.  Default is 1.
        min_inter_clade_dist : int
            The minimum inter-clade distance.  Default is 1.
        num_threads : int
            The number of threads desired for Ray processing. Default is 1.
        max_snp_resolution_thresh : int
            The maximum snp resolution for the clade.  The maximum number of members a clade is required to have to be considered valid.  Default is 0.
        method : str
            Method of clustering according to `scipy.cluster.hierarchy.linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage>`_.  Default is 'average'.
        rcor_thresh : float
            The correlation coefficient threshold.  Default is 0.4.
        min_perc : float
            Minimum percentage of clade members to be positive for a kmer to be valid.  Default is 1.
        """
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
        self.rcor_thresh = rcor_thresh
        self.max_snp_count = max_snps
        self.max_snp_resolution_thresh = max_snp_resolution_thresh
        self.min_perc = min_perc
        self.mode = mode
        self.workflow()

        return

    def workflow(self):
        """
        The workflow method calls all the helper methods of this class to process the data
        according to both the input and the flags set by the user to determine
        clade memberships, validate SNPs, and compress the resulting clade unless otherwise
        specified by the user.  Please refer to the file
        examples/small_test/cladeomatic/cladeomatic-clades.info.txt for more information
        """
        self.raw_genotypes = self.generate_genotypes()
        self.snp_data = snp_search_controller(self.group_data, self.vcf_file, self.num_threads)
        self.summarize_snps()
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
        #re-get the valid notes following the above processing
        valid_nodes = self.get_valid_nodes()
        self.set_valid_nodes(valid_nodes)
        self.dist_based_nomenclature()
        self.selected_genotypes = self.generate_genotypes()


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
            self.compression_cleanup()
        elif self.mode == 'tree':
            self.fix_root()
            valid_nodes = self.get_valid_nodes()
            self.set_valid_nodes(valid_nodes)

        for node_id in self.group_data['valid_nodes']:
            if node_id in self.clade_data:
                self.clade_data[node_id]['is_valid'] = True
        self.selected_genotypes = self.generate_genotypes()
        self.selected_positions = self.get_selected_positions()
        self.temporal_signal()
        self.set_genotype_snp_states()
        self.set_genotype_snp_rules()

    def update(self):
        """
        A method to call the helper methods to update the
        clades - :meth:`get_valid_nodes`, :meth:`set_valid_nodes`, :meth:`generate_genotypes`.
        """
        self.check_nodes()
        valid_nodes = self.get_valid_nodes()
        self.set_valid_nodes(valid_nodes)
        self.supported_genotypes = self.generate_genotypes()

    def summarize_snps(self):
        """
        A method to loop through the snps identified and determine
        if the snp is canonical and/or valid.  It also creates a list
        of variant positions for downstream processing.
        """
        num_canonical = 0
        num_valid = 0
        variant_pos = []
        num_positions = 0
        #loop through the snps identified from the vcf file
        for chrom in self.snp_data:
            num_positions = len(self.snp_data[chrom])
            for pos in self.snp_data[chrom]:
                variant_pos.append(pos)
                pos_is_valid = False
                pos_is_canonical = False
                for base in self.snp_data[chrom][pos]:
                    if self.snp_data[chrom][pos][base]['is_valid']:
                        pos_is_valid = True
                    if self.snp_data[chrom][pos][base]['is_canonical']:
                        pos_is_canonical = True
                if pos_is_valid:
                    num_valid+=1
                if pos_is_canonical:
                    num_canonical+=1
        self.num_positions = num_positions
        self.num_valid_positions = num_valid
        self.num_canonical_positions = num_canonical
        self.variant_positions = variant_pos

    def compression_cleanup(self):
        """
        A method to clean up the clade following compression by
        removing nodes whose distance is below the threshold default or
        the threshold as denoted by the user input.
        """
        node_below_threshold = set()
        #loop through the clade ids to find nodes below the threshold distance
        for clade_id in self.clade_data:
            dist = self.clade_data[clade_id]['ave_within_clade_dist']
            if dist < self.max_snp_resolution_thresh:
                node_below_threshold.add(clade_id)

        node_to_link_map = {}
        #remake the cluster node map, eliminating possible duplicates and optimizing the data structure
        for link_id in self.cluster_to_node_mapping:
            node_id = self.cluster_to_node_mapping[link_id]
            if not node_id in node_to_link_map:
                node_to_link_map[node_id] = set()
            node_to_link_map[node_id].add(link_id)

        #Remove nodes which are less than the distance threshold
        genotype_assignments = {}
        for sample_id in self.supported_genotypes:
            genotype = self.supported_genotypes[sample_id].split(self.delim)
            if len(genotype) <= 2:
                continue
            filt = []
            for node_id in genotype:
                if node_id in node_below_threshold:
                    continue
                filt.append(node_id)
            genotype_assignments[sample_id] = filt

        terminal_nodes = set()
        #determine which nodes are terminal nodes
        for sample_id in genotype_assignments:
            node_id = genotype_assignments[sample_id][-1]
            terminal_nodes.add(node_id)

        #remove the nodes below the threshold
        valid_nodes = self.get_valid_nodes() - node_below_threshold
        #add the terminal nodes as part of the valid nodes set
        valid_nodes = valid_nodes | terminal_nodes
        self.set_valid_nodes(valid_nodes)

    def get_selected_positions(self):
        """
        A helper method to retrieve the list of valid nodes within the clade data dictionary.

        Returns
        -------
        list
            A sorted list of integers for the positions of valid nodes
        """
        selected_positions = set()
        for node_id in self.clade_data:
            if self.clade_data[node_id]['is_valid']:
                selected_positions = selected_positions | set(self.clade_data[node_id]['pos'])
        return sorted(list(selected_positions))

    def set_valid_nodes(self,valid_nodes):
        """
        Set the valid nodes in the group data set.

        Parameters
        ----------
        set
            The valid nodes to set in the group data set
        """
        self.group_data['valid_nodes'] = set(valid_nodes)

    def set_invalid_nodes(self,invalid_nodes):
        """
        A helper method to loop through all the invalid nodes in the set
        and set these nodes as invalid in the clade data dictionary.

        Parameters
        ----------
        set
            The set of invalid nodes to mark as invalid
        """
        for node_id in invalid_nodes:
            self.clade_data[node_id]['is_valid'] = False
            self.clade_data[node_id]['is_selected'] = False

    def get_all_nodes(self):
        """
        Set the list of nodes as the valid nodes from the group data dictionary.
        """
        self.node_list = sorted(list(self.group_data['valid_nodes']))

    def get_valid_nodes(self):
        """
        A helper method to loop through all the nodes to find the
        nodes flagged as valid in the clade data set.

        Returns
        -------
        set
            A set of all the valid nodes in the clade data
        """
        valid_nodes = set()
        for node_id in self.clade_data:
            if self.clade_data[node_id]['is_valid']:
                valid_nodes.add(node_id)
        return valid_nodes

    def get_terminal_nodes(self):
        """
        A helper method to find the terminal nodes in the genotypes set generated in
        the :meth:`generate_genotypes` method.

        Returns
        -------
        set
            The set of terminal nodes from the genotypes set generated
        """
        genotypes = self.generate_genotypes()
        terminal_nodes = set()
        #retrieve the terminal node from the genotypes set
        for sample_id in genotypes:
            terminal_nodes.add(genotypes[sample_id].split(self.delim)[-1])

        return terminal_nodes

    def get_selected_nodes(self):
        """
        A helper method to retrieve the nodes that correspond to the 'is selected' flag within
        the clade data dictionary.

        Returns
        -------
        set
            The nodes that correspond to the 'is selected' flag in the clade data
        """
        valid_nodes = set()
        for node_id in self.clade_data:
            if self.clade_data[node_id]['is_selected']:
                valid_nodes.add(node_id)
        return valid_nodes

    def get_node_member_counts(self):
        """
        A helper method to count the number of distinct genotypes in the
        group data sample map and set that to the node counts class variable.
        """
        counts = {}
        for sample_id in self.group_data['sample_map']:
            genotype = self.group_data['sample_map'][sample_id]['genotype'].split(self.delim)
            for g in genotype:
                if g not in counts:
                    counts[g] = 0
                counts[g]+=1
        self.node_counts = counts

    def init_clade_data(self):
        """
        A method to initialize the clade data dictionary and all the
        variables contained within for each item in the node list data set.
        """
        for node_id in self.node_list:
            self.clade_data[node_id] = {
                'pos': [],
                'bases': [],
                'min_dist': 0,
                'max_dist': 0,
                'ave_dist': 0,
                'fisher': {},
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
                'clade_sample_id':None,
                'spearmanr':0,
                'spearmanr_pvalue':1,
                'pearsonr':0,
                'pearsonr_pvalue':1,
                'is_temporal_signal_present':False
            }

    def populate_clade_data(self):
        """
        The method for populating clade data dictionary with the
        canonical snps from the snp data set.  Uses the clade id,
        the position and variant base of the snp.
        """
        #loop through the bases of the snps identified
        for chrom in self.snp_data:
            for pos in self.snp_data[chrom]:
                for base in self.snp_data[chrom][pos]:
                    #if the snp is canonical, initialize the data in the clade data dictionary
                    if self.snp_data[chrom][pos][base]['is_canonical']:
                        clade_id = self.snp_data[chrom][pos][base]['clade_id'].split('.')[-1]
                        self.clade_data[clade_id]['pos'].append(pos)
                        self.clade_data[clade_id]['bases'].append(base)

    def check_nodes(self):
        """
        A method to check the validity of the various nodes in the clade data.
        If the number of nodes is below the :attr:`min_member_count` or
        the :attr:`min_snp_count`, the node flags 'is_valid' and 'is_selected' are
        set to False.
        """
        for node_id in self.clade_data:
            #the node must have greater than the minimum number of members to be considered valide
            if self.clade_data[node_id]['num_members'] < self.min_member_count:
                self.clade_data[node_id]['is_valid'] = False
                self.clade_data[node_id]['is_selected'] = False
            #the position must have greater than the minimum number of snps attached to be valid
            if len(self.clade_data[node_id]['pos']) < self.min_snp_count:
                self.clade_data[node_id]['is_valid'] = False
                self.clade_data[node_id]['is_selected'] = False

    def generate_genotypes(self):
        """
        A method to generate the genotypes data dictionary.  This method
        creates this set from the group data dictionary and extracts the
        tree node and genotype from the sample map contained in the group data dictionary.

        Returns
        -------
        dict
            A dictionary with the genotypes consisting of the tree node identifiers and genotype identifiers
        """
        #set a delimiter if none was defined
        if self.delim == None:
            self.delim = '.'

        genotypes = {}
        #loop through the sample map of the group data dictionary
        for id in self.group_data['sample_map']:
            sample_id = self.group_data['sample_map'][id]['sample_id']
            genotype = [str(x) for x in self.group_data['sample_map'][id]['genotype'].split(self.delim)]
            filt_genotype = []
            #filter for the valid nodes only
            for i in range(0, len(genotype)):
                node_id = genotype[i]
                if node_id in self.group_data['valid_nodes']:
                    filt_genotype.append(node_id)
            genotypes[sample_id] = "{}".format(self.delim).join(filt_genotype)
        return genotypes

    def snp_based_filter(self):
        """
        A helper method to filter the snps based on if they are valid
        based on the number of max states for the snp positions and the
        minimum member count.  If there are greater than the max states or fewer
        than the minimum member count, the snp is marked as invalid.
        """
        for chrom in self.snp_data:
            for pos in self.snp_data[chrom]:
                if len(self.snp_data[chrom][pos]) > self.max_states:
                    self.snp_data[chrom][pos][base]['is_valid'] = False
                for base in self.snp_data[chrom][pos]:
                    if self.snp_data[chrom][pos][base]['is_canonical']:
                        if self.snp_data[chrom][pos][base]['num_members'] < self.min_member_count:
                            self.snp_data[chrom][pos][base]['is_valid'] = False

    def get_bifurcating_nodes(self):
        """
        A method to split nodes with more than one genotype node identifier,
        the first is the parent and the second is the child node.  It then sets the
        results in the bifurcating node dictionary.
        """
        nodes = {}
        #unused
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
        """
        A helper method to retrieve the genotype from the sample id passed

        Parameters
        ----------
        sample_list : list
            A list of sample identifiers

        Returns
        -------
        dict
            The dictionary of the genotypes corresponding to the list of sample ids
        """
        lookup = {}
        genotypes = self.generate_genotypes()
        d = self.delim
        for sample_id in sample_list:
            lookup[sample_id] = genotypes[sample_id].split(d)

        return lookup

    def calc_node_distances(self):
        """
        This method reads the distance matrix file previously created and
        sets the node's closest distance, the closest sample distances, the average
        distances within the clade, total comparisons, and total distance
        for each node in the clade data dictionary.
        """
        fh = open(self.distance_matrix_file,'r')
        header = next(fh).rstrip().split("\t")
        num_columns = len(header)
        #get the sample id list from the header
        sample_list = header[1:]
        sample_lookup = self.genotype_lookup(sample_list)
        histo = {}
        index = 1
        #go through the distance matrix entries
        for line in fh:
            line = line.rstrip().split("\t")
            sample_id_1 = line[0]
            nodes_1 = set(sample_lookup[sample_id_1])
            #loop through the associated nodes
            for i in range(index,num_columns):
                sample_id_2 = header[i]
                #if the node id in the row is the same as the node id in the header, skip
                if sample_id_1 == sample_id_2:
                    continue
                #get the sample ids for the node
                nodes_2 = set(sample_lookup[sample_id_2])
                value = int(line[i])
                #set up the histogram frequency data
                if not value in histo:
                    histo[value] = 0
                histo[value] +=1
                #create a set with the values common to both nodes - overlap
                ovl = nodes_1 & nodes_2
                #set the distance for the node id and how many comparisons were made
                for node_id in ovl:
                    self.clade_data[node_id]['total_within_clade_dist'] += value
                    self.clade_data[node_id]['num_comparisons_within_clade_dist'] += 1
                #create a set with values that do not overlap with the second node
                no_ovl_1 =  nodes_1 - nodes_2
                for node_id in no_ovl_1:
                    if self.clade_data[node_id]['closest_sample_id'] is None or self.clade_data[node_id]['closest_sample_dist'] < value:
                        self.clade_data[node_id]['clade_sample_id'] = sample_id_1
                        self.clade_data[node_id]['closest_sample_id'] = sample_id_2
                        self.clade_data[node_id]['closest_sample_dist'] = value
                # create a set with values that do not overlap with the first node
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
            #calculate and set the average within clade distance of the samples
            if self.clade_data[node_id]['num_comparisons_within_clade_dist'] > 0:
                ave = self.clade_data[node_id]['total_within_clade_dist'] / self.clade_data[node_id]['num_comparisons_within_clade_dist']
                self.clade_data[node_id]['ave_within_clade_dist'] = ave
            sample_id_1 = self.clade_data[node_id]['clade_sample_id']
            sample_id_2 = self.clade_data[node_id]['closest_sample_id']
            #if there are no samples associated with the node identifier, skip
            if sample_id_1 is None:
                continue
            nodes_1 = sample_lookup[sample_id_1]
            nodes_2 = sample_lookup[sample_id_2]
            num_nodes = len(nodes_1)
            if len(nodes_2) < num_nodes:
                num_nodes = len(nodes_2)
            shared_parent = None
            #verify the rootings/parent nodes
            for i in range(0,num_nodes):
                n1 = nodes_1[i]
                n2 = nodes_2[i]
                if n1 == n2:
                    shared_parent = n1
            self.clade_data[node_id]['closest_clade_id'] = shared_parent
            self.clade_data[node_id]['closest_clade_dist'] = self.clade_data[node_id]['closest_sample_dist']

        self.distance_histo = histo

    def get_clade_distances(self):
        """
        This method loops though the clade data nodes for the average distances
        within the clades and places them in a list.

        Returns
        -------
        list
            A list of floats for the average distances within the clades
        """
        distances = []
        for node_id in self.clade_data:
            distances.append(self.clade_data[node_id]['ave_within_clade_dist'])
        return distances

    @deprecated
    def get_inter_clade_distances(self):
        """
        This helper method loops through the nodes of the clade data dictionary
        and creates a list of the closest clade distances for each node.

        Returns
        -------
        list
            A list of the closest clade distances
        """
        distances = []
        for node_id in self.clade_data:
            distances.append(self.clade_data[node_id]['closest_clade_dist'])
        return distances

    def get_close_nodes(self):
        """
        This helper method to find the nodes whose distances are smaller than
        the minimum inter-clade distance.

        Returns
        -------
        set
            A set for the invalid nodes who violate the distance constraint
        """
        invalid_nodes = set()
        for node_id in self.clade_data:
            dist = self.clade_data[node_id]['closest_clade_dist']
            if dist < self.min_inter_clade_dist:
                invalid_nodes.add(node_id)
                if self.clade_data[node_id]['closest_clade_id'] is not None:
                    invalid_nodes.add(self.clade_data[node_id]['closest_clade_id'])
        return invalid_nodes

    def partition_distances(self):
        """
        A helper method to determine the bins and partitions required for
        the clades based on the average distances between these clades.

        Returns
        -------
        list
            A list of the average distances (partitions) for the clade for the calculated bins
        """
        num_partitions = 10
        #get the average distances between the clades
        distances = sorted(self.get_clade_distances())
        num_distances = len(distances)
        #determine how many samples in a bin
        target_number_samples = int(num_distances / num_partitions)

        bins = {}
        #initialize the bins
        for i in range(0,num_partitions):
            bins[i] = []
        bin_index = 0
        #loop through the number of distances to add the distances to the
        #correct bins
        for i in range(0,num_distances):
            d = distances[i]

            if len(bins[bin_index]) > target_number_samples:
                if d > distances[i-1]:
                    bin_index+=1
            bins[bin_index].append(d)

        partitions = []
        #loop through the bins to add the last distance in the bin to the partition
        for i in bins:
            if len(bins[i]) > 0:
                partitions.append(bins[i][-1])

        return partitions

    def distance_node_binning(self):
        """
        A method to add the clade nodes to the bins as per the partition
        distances determined in the :meth:`partition_distances` method.
        """
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
        """
        A method to compress the generated hierarchy based on the node bins
        previously determined and the genotypes.  This method also removes nodes based
        on if the nodes on the same branch belong to the same bins.  Valid
        nodes are those nodes which are terminal and distinct.
        """
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
                    #if the nodes belong to the same bin, set the second node as invalid
                    if n1 in node_bins[bin_id] and n2 in node_bins[bin_id]:
                        invalid_nodes.add(n2)
                        break
        self.set_invalid_nodes(invalid_nodes)
        #remove the invalid nodes from the set
        self.set_valid_nodes(self.get_valid_nodes() - invalid_nodes)
        terminal_nodes = set()
        genotypes = self.generate_genotypes()
        #loop through to find the terminal nodes
        for sample_id in genotypes:
            genotype = genotypes[sample_id].split(self.delim)
            node_id = genotype[-1]
            if node_id in self.clade_data:
                terminal_nodes.add(node_id)
        #re-set the valid and invalid nodes based on the above compression
        invalid_nodes = invalid_nodes | (self.get_valid_nodes() - terminal_nodes)
        self.set_invalid_nodes(invalid_nodes)
        self.set_valid_nodes(terminal_nodes - invalid_nodes)

    def get_variant_positions(self):
        """
        A method to retrieve the variant positions from the snp data dictionary.

        Returns
        -------
        list
            An integer list of the variant positions
        """
        variant_positions = []
        for chrom in self.snp_data:
            for pos in self.snp_data[chrom]:
                variant_positions.append(pos)

        return variant_positions

    def get_conserved_ranges(self):
        """
        A method to compile a list of sequence ranges for the variant positions
        data list.

        Returns
        -------
        list
            A list of the conserved ranges for the variant positions
        """
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
        """
        A method to loop through the conserved sequence ranges determined through
        the method :meth:`get_conserved_ranges` and find the longest sequence and the
        position for the snp contained therein.

        Returns
        -------
        int
            The position for the snp within the largest conserved sequence
        """
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
        """
        A helper method to reset the root of the clade with the snp contained
        in largest conserved sequence.
        """
        pos = self.get_largest_conserved_pos()
        base = self.ref_seq[pos]
        self.clade_data['0']['pos'].append(pos)
        self.clade_data['0']['bases'].append(base)
        self.clade_data['0']['is_valid'] = True
        self.clade_data['0']['is_selected'] = True

    def perform_clustering(self):
        """
        A method to read the distance matrix file created, along with the
        distance thresholds, and perform the clustering though the
        SciPy hierarchy methods.

        Notes
        -----
        See `scipy.cluster.hierarchy.linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage>`_
        and `scipy.cluster.hierarchy.fcluster <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html>`_
        for more thorough documentation.
        """
        #read the distance matrix file and extract the values
        data = pd.read_csv(self.distance_matrix_file, sep='\t', header=0, index_col=0)
        data_labels = data.columns.values
        distance_matrix = distance.squareform(data.values)
        #create the hierarchy through the method denoted in the user specifications (default was average)
        Z = hierarchy.linkage(distance_matrix, method=self.method)
        cluster_membership = {}
        for id in data_labels:
            cluster_membership[id] = []
        clust_number = 1
        #loop through the distance thresholds to flatten the above clustering and cluster the samples
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
        """
        A method that loops through the sample map dictionary
        and should the sample id from the map match the cluster
        determined in the :meth:`perform_clustering` method, the node ids
        and sample map indexes are added to the nodes dictionary.

        Returns
        -------
        dict
            A dictionary with the sample map indexes and the node ids clustered
        """
        nodes = {}
        sample_map = self.group_data['sample_map']
        for index in sample_map:
            sample_id = sample_map[index]['sample_id']
            #retrieve the node from the clusters that corresponds
            #to the sample id from the group data and set that node
            #id and the sample index in the dictionary
            for node_id in self.sample_linkage_clusters[sample_id]:
                if not node_id in nodes:
                    nodes[node_id] = set()
                nodes[node_id].add(index)
        return nodes

    def fit_clusters_to_groupings(self):
        """
        A method to map the clusters created to the nodes previously identified.

        Returns
        -------
        dict
            A dictionary for the mapping of nodes to clusters
        """
        cluster_to_node_mapping = {}
        sl_cluster_membership = self.create_cluster_membership_lookup()
        #loop through the clusters in the cluster memberships
        for sl_clust_id in sl_cluster_membership:
            clust_members = sl_cluster_membership[sl_clust_id]
            best_node_id = '0'
            best_node_nonmembers = len(self.group_data['membership'])
            #loop through the group data membership nodes
            for node_id in self.group_data['membership']:
                node_members = self.group_data['membership'][node_id]
                #find the node members unique to the cluster memberships
                missing = clust_members - node_members
                #skip the next part if there are missing node members in the cluster
                if len(missing) > 0:
                    continue
                #find the nonmembers that are unique to the group data membership
                nonmembers = node_members - clust_members
                count_nonmembers = len(nonmembers)
                #if the count of the non-members is smaller than the best count
                #set the variables to find the group that has the least difference
                #between sets
                if count_nonmembers < best_node_nonmembers:
                    best_node_id = node_id
                    best_node_nonmembers = count_nonmembers
                if count_nonmembers == 0:
                    break
            #map the cluster id to the best node members
            cluster_to_node_mapping[sl_clust_id] = best_node_id.split(self.delim)[-1]
        return cluster_to_node_mapping

    def identify_dist_ranges(self):
        """
        A method to identify the troughs within the histogram data.  This method
        uses the helper method :meth:`find_dist_troughs` to accomplish this task.

        Returns
        -------
        list
            A ist of tuples for the range of the troughs in the histogram data
        """
        #order the histogram dictionary
        histo = OrderedDict(sorted(self.distance_histo.items()))
        #retrieve the keys
        histo_keys = list(histo.keys())
        max_value = histo_keys[-1]
        x = list(range(0,max_value+2))
        y = []
        #loop through and get the values for the keys
        for i in x:
            value = 0
            if i in histo:
                value = histo[i]
            y.append(value)
        #find the troughs
        troughs, properties = self.find_dist_troughs(y)
        self.dist_thresholds = sorted(list(troughs),reverse=True)
        troughs = sorted(list(set(troughs)))
        c = self.as_range(troughs)
        ranges = []
        #set the ranges
        for idx, tuple in enumerate(c):
            ranges.append((troughs[tuple[0]], troughs[tuple[1]]))

        return ranges

    def dist_based_nomenclature(self):
        """
        A method to find the genotype identifiers from clustered data
        and format these identifiers based on the distances found and the delimiter set
        upon :class:`clade_worker` instantiation.
        """
        #perform the clustering
        self.dist_ranges = self.identify_dist_ranges()
        self.perform_clustering()
        clust_to_node_mapping = self.fit_clusters_to_groupings()
        self.cluster_to_node_mapping = clust_to_node_mapping
        sample_genotypes = {}
        for sample_id in self.sample_linkage_clusters:
            clust_geno = self.sample_linkage_clusters[sample_id]
            node_geno = []
            #map the cluster to the nodes
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
        self.sample_dist_based_nomenclature = sample_genotypes

    def find_dist_troughs(self,values):
        """
        Uses the SciPy method `scipy.signal.find_peaks <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html>`_
        to locate the troughs (reverse peaks) in the histogram data passed.

        Parameters
        ----------
        values : list
            The list of the histogram data values to find the troughs for

        Returns
        -------
        troughs : numpy.array
            An array of the troughs
        properties :dict
            a dictionary of the properties

        Notes
        _____
        Please refer to the external methods `scipy.signal.find_peaks <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html>`_
        and `numpy.array <https://numpy.org/doc/stable/reference/generated/numpy.array.html>`_
        for more thorough documentation.

        """
        return find_peaks(-np.array(values))

    def as_range(self,values):
        """
        A helper method to determine the trough ranges through the histogram data passed.

        Parameters
        ----------
        values : list
            A list of integers to find the ranges for

        Returns
        -------
        list
            A list of tuples with the start and end index of each range
        """
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
        """
        This function removes snps from the valid positions variable in
        this clade data dictionary.

        Parameters
        ----------
        positions : list
            A list of integer positions for the snps to be removed

        """
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

    def calc_metadata_counts(self,sample_set):
        """
        A helper method to count the metadata values in the sample set passed.
        This method loops through the field identifiers of the
        metadata and counts the total values for each field value.

        Example: the location field has the values Europe and North America.  This
        method will count the number of samples that correspond to Europe and North
        America (Europe 10, North America 15).

        Parameters
        ----------
        sample_set : set
            A set for the samples to tabulate metadata counts

        Returns
        -------
        dict
            A dictionary for the metadata field id, value and the value counts
        """
        metadata_counts = {}
        #find the sample ids in the metadata dictionary that correspond to the sample set
        for sample_id in self.metadata_dict:
            if sample_id not in sample_set:
                continue
            #loop through the field identifiers in the metadata set and count
            #the number of times the value appears
            for field_id in self.metadata_dict[sample_id]:
                value = self.metadata_dict[sample_id][field_id]
                if not field_id in metadata_counts:
                    metadata_counts[field_id] = {}
                if not value in metadata_counts[field_id]:
                    metadata_counts[field_id][value] = 0
                metadata_counts[field_id][value] +=1
        return metadata_counts

    def calc_node_associations_groups(self):
        """
        A method to loop through the metadata counts to
        calculate the fisher's exact test for the node associations
        in the clades aleady determined.
        """
        #retrieve the sample ids
        samples = set(self.metadata_dict.keys())
        num_samples = len(samples)
        #count the field values in the sample metadata
        metadata_counts = self.calc_metadata_counts(samples)
        #loop through and initialize each of the field names for the fisher exact test
        for clade_id in self.clade_data:
            for col in metadata_counts:
                if col == 'year':
                    continue
                for field_name in metadata_counts[col]:
                    self.clade_data[clade_id]['fisher'][field_name] = {'oddsr': 'nan', 'p': 1}

        #loop through the genotypes of the group member data to get the membership and metadata field counts
        for genotype in self.group_data['membership']:
            clade_id = genotype.split(self.delim)[-1]
            in_members_indexes = self.group_data['membership'][genotype]
            in_members = set()
            for id in in_members_indexes:
                in_members.add(self.group_data['sample_map'][id]['sample_id'])
            num_members = len(in_members)
            in_counts = self.calc_metadata_counts(in_members)
            for col in in_counts:
                if col == 'year':
                    continue
                #construct a pseudo table with the metadata counts, the number of members
                #and the number of samples to calculate the fisher's exact test coefficient
                #between fields - to determine the association between them, if any
                for field_name in in_counts[col]:
                    table  = [
                        [ 0, 0 ],
                        [ 0, 0 ],
                    ]
                    pp = in_counts[col][field_name]
                    pn = num_members - pp
                    np = metadata_counts[col][field_name] - pp
                    nn = num_samples - metadata_counts[col][field_name]  - pn
                    #heuristic to skip results with low likelihood of success
                    if pn / num_members > 0.2:
                        continue
                    table[0][0] = pp
                    table[0][1] = pn
                    table[1][0] = np
                    table[1][1] = nn
                    oddsr, p = fisher_exact(table, alternative='greater')
                    self.clade_data[clade_id]['fisher'][field_name] = {'oddsr': oddsr, 'p': p}

    def temporal_signal(self):
        """
        A method to calculate the spearman's and pearson's coefficients
        from the metadata for the year and the calculated clade distances
        previously determined.  If the spearman or pearson coefficients are
        larger than the pre-determined threshold, the temporal flag is set
        to 'True' in the clade data for the node.
        """
        #find the year of the sample in the metadata if it exists
        sample_id = list(self.metadata_dict.keys())[0]
        if 'year' not in self.metadata_dict[sample_id]:
            return

        temporal_data = {}
        #initialize the temporal data dictionary
        for node_id in self.clade_data:
            temporal_data[node_id] = {
                'dist':[],
                'year':[]
            }
        #open the matrix file
        fh = open(self.distance_matrix_file, 'r')
        header = next(fh).rstrip().split("\t")
        num_columns = len(header)
        sample_list = header[1:]
        #get the genotypes for the samples
        sample_lookup = self.genotype_lookup(sample_list)
        index = 1
        #loop through the distance matrix rows
        for line in fh:
            line = line.rstrip().split("\t")
            #get the sample ids and member nodes
            sample_id_1 = line[0]
            nodes_1 = set(sample_lookup[sample_id_1])
            y1 = self.metadata_dict[sample_id_1]['year']
            if y1 == 'nan' or len(y1) == 0 or y1 is None:
                continue
            y1 = float(y1)
            #loop through to find the distance between two samples and ensure they both have years
            #set the temporal data with the node id as key and the distance and year value
            for i in range(index, num_columns):
                sample_id_2 = header[i]
                if sample_id_1 == sample_id_2:
                    continue
                y2 = self.metadata_dict[sample_id_2]['year']
                if y2 == 'nan' or len(y2) == 0 or y2 is None:
                    continue
                y2 = float(y2)
                nodes_2 = set(sample_lookup[sample_id_2])
                value = int(line[i])
                ovl = nodes_1 & nodes_2
                for node_id in ovl:
                    temporal_data[node_id]['dist'].append(value)
                    temporal_data[node_id]['year'].append(y2)
            for node_id in ovl:
                temporal_data[node_id]['dist'].append(value)
                temporal_data[node_id]['year'].append(y1)
        fh.close()

        #loop through the temporal data created above
        for node_id in temporal_data:
            R = 0
            P = 1
            Rpear = 0
            Ppear = 1
            num_distinct_dist = len(set(temporal_data[node_id]['dist']))
            num_distinct_year = len(set(temporal_data[node_id]['year']))
            #calculate the spearman's and pearson's coefficients for the temporal
            #data based on the distance and year
            if num_distinct_dist >= 3 and num_distinct_year >= 3:
                R, P = spearmanr(np.asarray(temporal_data[node_id]['year']), np.asarray(temporal_data[node_id]['dist']))
                Rpear, Ppear = pearsonr(np.asarray(temporal_data[node_id]['year']), np.asarray(temporal_data[node_id]['dist']))

            #set these values in the clade data for the node
            self.clade_data[node_id]['spearmanr'] = R
            self.clade_data[node_id]['spearmanr_pvalue'] = P
            self.clade_data[node_id]['pearsonr'] = Rpear
            self.clade_data[node_id]['pearsonr_pvalue'] = Ppear
            is_present = False
            #flag if there is temporal signal in the data based on the
            #spearman's and pearson's coefficients calculated above
            if R > self.rcor_thresh or Rpear > self.rcor_thresh:
                is_present = True
            self.clade_data[node_id]['is_temporal_signal_present'] = is_present

    def clade_snp_count(self):
        """
        A method to aggregate the counts for the number of clade node
        members based on the clade id in the clade data dictionary.

        Returns
        -------
        dict
            A dictionary of the counts for the clade node members
        """
        counts = {}
        for node_id in self.clade_data:
            counts[node_id] = len(self.clade_data[node_id]['pos'])
        return counts

    def clade_snp_association(self):
        """
        A helper method to loop through the clade and find the associated
        snp names for each clade.

        Returns
        -------
        dict
            A dictionary for the snp name and associated clade node id
        """
        snps = {}
        for node_id in self.clade_data:
            for idx,pos in enumerate(self.clade_data[node_id]['pos']):
                base = self.clade_data[node_id]['bases'][idx]
                name = "{}-{}".format(pos,base)
                if not name in snps:
                    snps[name] = set()
                snps[name].add(node_id)
        return snps

    def prune_snps(self):
        """
        A method to prune the nodes in the clade based on the :attr:`max_snps` constraint.
        If there are more snp entries for the clade node, then some are removed
        randomly.
        """
        max_snps = self.max_snp_count
        node_counts = self.clade_snp_count()
        nodes_to_prune = []

        for node_id in node_counts:
            count = node_counts[node_id]
            if count > max_snps:
                nodes_to_prune.append(node_id)

        snp_assoc = self.clade_snp_association()
        selected_snps = set()
        for name in snp_assoc:
            (pos,base) = name.split('-')
            pos = int(pos)
            if len(snp_assoc[name]) == 1:
                selected_snps.add(pos)

        for node_id in self.clade_data:
            if node_counts[node_id] < max_snps:
                continue
            #determine if there is overlap between the clade-snp association
            #and the selected snps
            no_ovl = set(self.clade_data[node_id]['pos']) - selected_snps
            if len(no_ovl) == 0:
                continue
            ovl = set(self.clade_data[node_id]['pos'])  - no_ovl
            pos = []
            bases = []
            #if the overlapping snps list is larger than the max snps constraint
            #shuffle the set and remove the required amount of snps for the
            #constraint
            if len(ovl) < max_snps:
                c = list(no_ovl)
                random.shuffle(c)
                ovl = ovl |set(c[0:max_snps-len(ovl)])
            for idx, p in self.clade_data[node_id]['pos']:
                if p in ovl:
                    pos.append(p)
                    bases.append(self.clade_data[node_id]['bases'][p])

            self.clade_data[node_id]['pos'] = pos
            self.clade_data[node_id]['bases'] = bases

        self.selected_positions = self.get_selected_positions()

    def set_genotype_snp_states(self):
        """
        This method sets the genotype data, genotype and base counts,
        the snp position, and in the genotype data for the
        :attr:`genotype_snp_data` data dictionary.
        """
        selected_positions = self.selected_positions
        vcf = vcfReader(self.vcf_file)
        data = vcf.process_row()
        samples = vcf.samples
        sample_genotypes = self.selected_genotypes
        genotype_data = {}
        for sample_id in sample_genotypes:
            genotype = sample_genotypes[sample_id]
            if not genotype in genotype_data:
                genotype_data[genotype] = {
                    'genotype_count':0,
                    'base_counts':{}
                }
            genotype_data[genotype]['genotype_count']+=1
            for pos in selected_positions:
                genotype_data[genotype]['base_counts'][pos] = {
                    'A': 0,
                    'T': 0,
                    'C': 0,
                    'G': 0,
                    'N': 0
                }

        if data is None:
            return {}

        while data is not None:
            ref = data['REF']
            pos = int(data['POS']) - 1
            if pos not in selected_positions:
                data = vcf.process_row()
                continue
            self.selected_pos_ref_states[pos] = ref
            for sample_id in samples:
                base = data[sample_id]
                genotype = sample_genotypes[sample_id]
                if base not in ['A', 'T', 'C', 'G']:
                    base = 'N'
                genotype_data[genotype]['base_counts'][pos][base] += 1
            data = vcf.process_row()

        self.genotype_snp_data = genotype_data
        del (vcf)


    def set_genotype_snp_rules(self):
        """
        This method sets the genotype and snp rules for what is a
        positive and partial genotype within the genotype_snp associated
        data dictionary based on what the minimum percentage of clade members
        that need to be positive for a kmer to be valid.
        """
        valid_bases = ["A","T","C","G"]
        rules = {}
        for pos in self.selected_positions:
            rules[pos] = {}
            for base in valid_bases:
                rules[pos][base] = {'positive_genotypes': [], 'partial_genotypes': []}

            for genotype in self.genotype_snp_data:
                total = self.genotype_snp_data[genotype]['genotype_count']
                if total == 0:
                    continue
                for base in valid_bases:
                    p = self.genotype_snp_data[genotype]['base_counts'][pos][base] / total
                    if p >= self.min_perc:
                        rules[pos][base]['positive_genotypes'].append(genotype)
                    elif p > 0:
                        rules[pos][base]['partial_genotypes'].append(genotype)
        self.genotype_snp_rules = rules

    def get_genotype_snp_rules(self):
        """
        Getter method to retrieve the :attr:`genotype_snp_rules` data dictionary.

        Returns
        -------
        dicts
            The :attr:`genotype_snp_rules` data dictionary
        """
        return self.genotype_snp_rules
