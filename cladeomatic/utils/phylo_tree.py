import os

import dendropy
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, RectFace
from deprecated import deprecated


if not 'DISPLAY' in os.environ:
    os.environ['QT_QPA_PLATFORM']='offscreen'

@deprecated()
def tree_to_distance_matrix(tree_file,out_file):
    """
    This method converts a dendropy tree object to a distance
    matrix and writes this matrix to a csv file.

    Parameters
    ----------
    tree_file : str
        The path to the tree file (newick format only)
    out_file : str
        The path to the CSV output file for the distance matrix
    """
    tree = dendropy.Tree.get(path=tree_file, schema='newick')
    pdm = tree.phylogenetic_distance_matrix()
    pdm.write_csv(out_file)

@deprecated()
def get_pairwise_distances_from_matrix(matrix_file):
    """
    This method gets the pairwise distances from the passed matrix file
    and adds them to a list.

    Parameters
    ----------
    matrix_file : str
        The file path to the distance matrix file
    Returns
    -------
    list
      The list of the pairwise distances
    """
    dists = []
    with open(matrix_file,'r') as fh:
        line = next(fh)
        for line in fh:
            row = line.split(',')
            dists += [float(x) for x in row[1:]]
        fh.close()
    return dists

def remove_unsupported_clades(ete_tree_obj, valid_clades):
    """
    Remove non-leaf nodes from the tree while maintaining the rest of the child nodes in the tree.

    Parameters
    ----------
    ete_tree_obj:
        The phylogenetic tree created by Clade-o-matic
    valid_clades : list
        The list of the valid nodes to maintain in the tree

    Returns
    -------
    ete_tree_obj: ETE3 tree object
        The phylogenetic tree with the invalid nodes removed

    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation.
    """
    tree = ete_tree_obj.copy()

    for node in tree.traverse("preorder"):
        if node.is_root() or node.is_leaf():
            continue

        node_id = node.name
        remove_node = True
        if node_id in valid_clades:
            remove_node = False
        if remove_node:
            node.delete()
    return tree

@deprecated()
def prune_tree(ete_tree_obj,valid_nodes):
    """
    This method removes nodes from the ETE3 tree object
    that are not part of the list of valid node names.

    Parameters
    ----------
    ete_tree_obj : ETE3 tree object
        The ETE3 tree object for pruning
    valid_nodes : list
        The list of node names which are valid and not to be removed

    Returns
    -------
    list
        A list of both the ETE3 tree object and the list of only valid nodes
    
    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation.
    """
    invalid_nodes = []
    #traverse the tree to make a list of invalid nodes
    for node in ete_tree_obj.traverse("postorder"):
        node_id = node.name
        if node.is_leaf() or node_id not in valid_nodes:
            continue
        node_ancestors = node.get_ancestors()

        for n in node_ancestors:
            ancestor_id = n.name
            if n.is_leaf() or \
                    ancestor_id not in valid_nodes or \
                    ancestor_id in invalid_nodes:
                continue
            children = n.get_children()
            count_internal_nodes = 0
            for c in children:
                if c.is_leaf() or c.is_root() or c.name in invalid_nodes:
                    continue
                count_internal_nodes += 1
            if count_internal_nodes < 1:
                invalid_nodes.append(c.name)
    valid_nodes = list(set(valid_nodes) - set(invalid_nodes))
    pruned_tree = remove_unsupported_clades(ete_tree_obj, valid_nodes)
    return [pruned_tree,valid_nodes]


def parse_tree(tree_file,logging,ete_format=0,set_root=False,resolve_polytomy=True,ladderize=True,method='midpoint',outgroup=''):
    """
    This method parses the newick formatted tree and transforms it to an ETE3 tree object.

    Parameters
    ----------
    tree_file :str
        The path to newick formatted tree for transformation
    logging : logging object
        The logging object for log output
    ete_format : int
        The integer code for the format of the ETE3 tree object.  Default is 0.  Please refer to the ETE3 documentation http://etetoolkit.org/docs/latest/reference/reference_tree.html
    set_root : bool
        A boolean to flag if the root is set for the tree. Default is False.
    resolve_polytomy :bool
        Force bifurcations of the tree on polytomies.  Default is True.
    ladderize : bool
        Sort the branches of a given tree (swapping children nodes) according to the size.  Default is True.
    method : str
        Method to root the resulting tree, either midpoint or outgroup.  Default is midpoint.
    outgroup : str
        Name of taxon to root tree on.  Default is a blank string.

    Returns
    -------
    ETE3 Tree Object
        The ETE3 tree object filled with the input data

    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation.
    """

    logging.info("Attempting to parse tree file {}".format(tree_file))

    if not os.path.isfile(tree_file):
        logging.error("Specified tree file {} is not found, please check that it exists".format(tree_file))
        return dict()
    if os.path.getsize(tree_file) == 0:
        logging.error("Specified tree file {} is found but is empty".format(tree_file))
        return dict()


    # Load a tree structure from a newick file.
    t = Tree(tree_file,format=ete_format)

    logging.info(
        "Read {} samples comprising {} nodes from tree {}".format(len(t.get_leaf_names()),
                                                                  len(t.get_descendants()), tree_file))

    #Need this otherwise groups derived from the tree are inaccurate
    if resolve_polytomy:
        logging.info("Resolving polytomies through forced bifurcation {}".format(tree_file))
        t.resolve_polytomy()


    #Ladderize tree for aesthetics
    if ladderize:
        logging.info("Ladderizing tree {}".format(tree_file))
        t.ladderize()


    #Rooting tree based on user criteria
    if set_root:
        if method == 'midpoint':
            logging.info("Setting root based on midpoint rooting from ETE3")
            root = t.get_midpoint_outgroup()
            t.set_outgroup(root)
        elif method == 'outgroup':
            if outgroup == '':
                logging.error("User selected outgroup rooting but did not provide an outgroup, please refer to the documentation for setting an outgroup")
                return None
            else:
                #if outgroup label not in the tree, return an error to the user
                if t.search_nodes(name=outgroup):
                    t.set_outgroup(outgroup)
                else:
                    logging.error(
                        "Specified outgroup not found in tree, please check the name and try again")
                    return None

    sample_list = t.get_leaf_names()

    #label all internal nodes
    node_id = 0
    for node in t.traverse("preorder"):
        if node.name == '':
            while node_id in sample_list:
                node_id+=1
            node.name = str(node_id)
            node_id+=1


    num_samples = len(t.children[0].get_tree_root())
    num_nodes = len(t.get_descendants())
    is_rooted = len(t.get_tree_root().children) == 2
    root_id = t.get_tree_root().name
    max_node, max_dist = t.get_farthest_leaf()
    logging.info(
        "Read {} samples comprising {} nodes from tree {}:\n".format(num_samples,num_nodes, tree_file,sample_list))
    logging.info("Most distant sample id {}".format(max_node.name))
    logging.info("Max Distance {}".format(max_dist))
    if is_rooted:
        for node in t.traverse("levelorder"):
            if node.is_leaf():
                root_id = node.name
                break
        logging.info("Tree is rooted on sample {}".format(root_id))
    else:
        logging.error("Error the tree is unrooted, you must either specify the root using an outgroup or midpoint rooting")

    return t


def get_internal_clades(ete_tree_obj):
    """
    Identify internal nodes within the ETE3 tree and return a dictionary of the samples
    associated with that tree node's clade.

    Parameters
    ----------
    ete_tree_obj : ETE3 tree object
        The ETE3 tree object for parsing

    Returns
    -------
    clade_dict : dict
        The dictionary of tree nodes which are not leaves

    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation.
    """
    cache = ete_tree_obj.get_cached_content()
    clade_dict = {}
    for clade in cache:
        clade_id = clade.name
        children = clade.get_leaf_names()
        if len(children) == 1:
            continue
        clade_dict[clade_id] = children
    return clade_dict

@deprecated()
def init_clade_info(ete_tree_obj):
    """
    This method initializes the clade dictionary for the tree passed.
    It parses the tree to find the internal clade memberships and
    sets these to a dictionary.

    Parameters
    ----------
    :param ete_tree_obj: ETE3 tree object
        The tree to parse for clade memberships

    Returns
    -------
    dict
        The initialized clade membership dictionary for the tree passed

    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation.
    """
    clade_info = {}
    memberships = get_internal_clades(ete_tree_obj)

    for clade_id in memberships:
        clade_info[clade_id] = {'num_members':len(memberships[clade_id]),'membership':set(memberships[clade_id]),'canSNPs':{},'canSNP_count':0}

    return clade_info


@deprecated()
def get_tree_node_distances(ete_tree_obj):
    """
    This method creates a dictionary of tree node distances for the
    ETE3 tree object passed.

    Parameters
    ----------
    ete_tree_obj: ETE3 tree object
        The tree to determine distances between nodes for

    Returns
    --------
    :return: dictionary - the dictionary of the node name and the distance calculated

    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation.
    """
    distances = {}
    for node in ete_tree_obj.iter_descendants("preorder"):
        distances[node.name] = node.dist
    return distances

@deprecated()
def get_tree_node_bootstrap(ete_tree_obj):
    """
    A method to determine the branch support for the branches in the
    ETE3 tree object passed.

    Parameters
    ----------
    ete_tree_obj : ETE3 tree object
        The tree to parse

    Returns
    -------
    dict
        The dictionary of node names and their support value

    Notes
    -----
    Refer to http://etetoolkit.org for more thorough ETE3 documentation.
    """
    support = {}
    for node in ete_tree_obj.iter_descendants("preorder"):
        support[node.name] = node.support

    return support

