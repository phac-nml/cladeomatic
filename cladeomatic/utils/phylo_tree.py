import os

import dendropy
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, RectFace


if not 'DISPLAY' in os.environ:
    os.environ['QT_QPA_PLATFORM']='offscreen'

def tree_to_distance_matrix(tree_file,out_file):
    '''

    Parameters
    ----------
    tree_file
    out_file

    Returns
    -------

    '''
    tree = dendropy.Tree.get(path=tree_file, schema='newick')
    pdm = tree.phylogenetic_distance_matrix()
    pdm.write_csv(out_file)

def get_pairwise_distances_from_matrix(matrix_file):
    '''

    Parameters
    ----------
    matrix_file: str path to file

    Returns list of distances
    -------

    '''
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
    Remove non-leaf nodes from the tree while maintaining children
    :param ete_tree_obj: [ETE obj] Phylogenetic tree
    :param valid_clades: [list] Valid node names to maintain
    :return: ete_tree_obj: [ETE obj] Phylogenetic tree with only nodes in the valid list
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

def prune_tree(ete_tree_obj,valid_nodes):
    """
    :param ete_tree_obj: [ETE obj] Tree object
    :param valid_nodes: [List] Node names which are valide
    :return: [ETE Tree ] With only valid nodes
    """
    invalid_nodes = []
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
    Parses newick formatted tree into an ETE tree object
    Parameters
    ----------
    tree_file [str] : Path to newick formatted tree
    logging [logging obj] : logging object
    ete_format [int] : Refer to documentation http://etetoolkit.org/docs/latest/reference/reference_tree.html
    set_root [Bool] : Set root for the tree
    resolve_polytomy [Bool]: Force bifurcations of the tree on polytomies
    ladderize [Bool] : Sort the branches of a given tree (swapping children nodes) according to the size
    method [str]: Method to root tree, either midpoint or outgroup
    outgroup [str] : Name of taxon to root tree on

    Returns
    -------

    ETE tree obj

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
    Identify internal nodes in the ETE3 tree and return a dict of the samples associated with that clade
    Parameters
    ----------
    ete_tree_obj [ETE tree obj] :

    Returns
    -------
    clade_dict [dict]: Dictionary of tree nodes which are not leaves

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

def init_clade_info(ete_tree_obj):
    clade_info = {}
    memberships = get_internal_clades(ete_tree_obj)

    for clade_id in memberships:
        clade_info[clade_id] = {'num_members':len(memberships[clade_id]),'membership':set(memberships[clade_id]),'canSNPs':{},'canSNP_count':0}

    return clade_info


def get_tree_node_distances(ete_tree_obj):
    distances = {}
    for node in ete_tree_obj.iter_descendants("preorder"):
        distances[node.name] = node.dist
    return distances

def get_tree_node_bootstrap(ete_tree_obj):
    support = {}
    for node in ete_tree_obj.iter_descendants("preorder"):
        support[node.name] = node.support

    return support

