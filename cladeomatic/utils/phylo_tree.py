import os

import dendropy
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, RectFace

from cladeomatic.constants import COLORS

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

def annotate_tree(outfile,ete_tree_obj,node_ids,leaf_meta={},node_colors=None,circular=True,font_size=18,show_branch_support=False,h=1600,w=900,dpi=1000):

    #fix issue with rooting causing rendering of circular trees to fail
    root_node = ete_tree_obj.get_tree_root()
    if root_node.dist == 0:
        root_node.dist = 0.00000000000000001

    node_ids = sorted(list(node_ids))
    #Autoscale content
    num_samples = len(ete_tree_obj.get_leaves())
    num_fields = 0
    field_lens = {}
    for id in leaf_meta:
        n = len(leaf_meta[id])
        if n > num_fields:
            num_fields = n
        fields = leaf_meta[id]

        for i in range(0,len(fields)):
            fields[i] = str(fields[i])
            l = len(fields[i])
            if not i in field_lens:
                field_lens[i] = 0
            if l > field_lens[i]:
                field_lens[i] = l
    h = num_samples * 10
    w = h

    #Set color pallet
    color_list = []
    max_secondary = 0
    for primary_color in COLORS:
        num = len(COLORS[primary_color])
        if num > max_secondary:
            max_secondary = num
    for i in range(0, max_secondary):
        for primary_color in COLORS:
            num = len(COLORS[primary_color])
            if i < num:
                secondary_color = list(COLORS[primary_color].keys())[i]
                color_list.append(COLORS[primary_color][secondary_color])

    # Basic tree style
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = show_branch_support
    ts.draw_aligned_faces_as_table = True
    if circular:
        ts.mode='c'


    # Init node backgrounds
    num_colors = len(color_list)
    node_colors = {}
    num_ranks = 0
    k = 0
    p = 0
    for sample_id in leaf_meta:
        genotype = leaf_meta[sample_id][0].split('.')
        l = len(genotype)
        if l > num_ranks:
            num_ranks = l
        for node in genotype:
            if node in node_colors:
                continue
            if k >= num_colors:
                k = 0
                p += 1
            node_colors[node] = color_list[k]
            k += 1

    # Applies the same static style to all nodes in the tree. Note that,
    # if "nstyle" is modified, changes will affect to all nodes
    for n in ete_tree_obj.traverse():
        name = n.name
        nstyle = NodeStyle()
        nstyle["vt_line_color"] = "black"
        nstyle["hz_line_color"] = "black"
        nstyle["vt_line_type"] = 0
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["size"] = 0
        nstyle["shape"] = "square"
        #if name in node_ids:
        #    nstyle["bgcolor"] = node_colors[name]

        n.set_style(nstyle)

    for n in ete_tree_obj.get_leaves():
        name = n.name

        index=1
        genotype = leaf_meta[name][0].split('.')
        glen = len(genotype) - 1
        label = leaf_meta[name][0]
        face = TextFace(ftype='Courier', text=label, fsize=font_size)
        face.hz_align = 0
        face.background.color = 'white'
        face.margin_right = 0
        face.margin_left = 1
        face.border.width = 1
        face.border.color = 'white'
        n.add_face(face, column=0, position="aligned")

        for i in range(0, num_ranks+1):

            if i > glen:
                face = RectFace(width=50, height=50, bgcolor='white',fgcolor='white')
                face.border.color = 'white'
            else:
                face = RectFace(width=50, height=50, bgcolor=node_colors[genotype[i]],fgcolor=node_colors[genotype[i]])
                face.border.color = node_colors[genotype[i]]
            face.margin_right = 0
            face.margin_left = 0
            face.border.width = 1


            face.inner_border.width = 0
            n.add_face(face, column=index , position="aligned")
            index+=1


    ete_tree_obj.render(outfile,tree_style=ts,dpi=dpi)

def plot_single_rep_tree(outfile,ete_tree_obj,node_ids,leaf_meta={},node_colors=None,circular=False,font_size=18,show_branch_support=False,h=1600,w=900,dpi=1000):
    #fix issue with rooting causing rendering of circular trees to fail
    root_node = ete_tree_obj.get_tree_root()
    if root_node.dist == 0:
        root_node.dist = 0.00000000000000001

    # Set color pallet
    color_list = []
    max_secondary = 0
    for primary_color in COLORS:
        num = len(COLORS[primary_color])
        if num > max_secondary:
            max_secondary = num
    for i in range(0, max_secondary):
        for primary_color in COLORS:
            num = len(COLORS[primary_color])
            if i < num:
                secondary_color = list(COLORS[primary_color].keys())[i]
                color_list.append(COLORS[primary_color][secondary_color])

    # Basic tree style
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = show_branch_support
    ts.margin_left = 10
    ts.margin_right = 10
    ts.margin_top = 10
    ts.margin_bottom = 10
    ts.branch_vertical_margin = 1
    ts.children_faces_on_top = False
    if circular:
        ts.mode='c'

    # Init node backgrounds
    num_colors = len(color_list)
    node_colors = {}
    num_ranks = 0
    k = 0
    p = 0
    for sample_id in leaf_meta:
        genotype = leaf_meta[sample_id][0].split('.')
        l = len(genotype)
        if l > num_ranks:
            num_ranks = l
        for node in genotype:
            if node in node_colors:
                continue
            if k >= num_colors:
                k = 0
                p += 1
            node_colors[node] = color_list[k]
            k += 1

    # Applies the same static style to all nodes in the tree. Note that,
    # if "nstyle" is modified, changes will affect to all nodes
    for n in ete_tree_obj.traverse():
        name = n.name
        nstyle = NodeStyle()
        nstyle["vt_line_color"] = "black"
        nstyle["hz_line_color"] = "black"
        nstyle["vt_line_type"] = 0
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["size"] = 0
        nstyle["shape"] = "square"
        # if name in node_ids:
        #    nstyle["bgcolor"] = node_colors[name]

        n.set_style(nstyle)

    for n in ete_tree_obj.get_leaves():
        name = n.name

        index = 1
        genotype = leaf_meta[name][0].split('.')
        glen = len(genotype) - 1
        label = leaf_meta[name][0]
        face = TextFace(ftype='Courier', text=label, fsize=font_size)
        face.hz_align = 0
        face.background.color = 'white'
        face.margin_right = 0
        face.margin_left = 1
        face.border.width = 1
        face.border.color = 'white'
        n.add_face(face, column=0, position="aligned")

        for i in range(0, num_ranks + 1):

            if i > glen:
                face = RectFace(width=50, height=50, bgcolor='white', fgcolor='white')
                face.border.color = 'white'
            else:
                face = RectFace(width=50, height=50, bgcolor=node_colors[genotype[i]],
                                fgcolor=node_colors[genotype[i]])
                face.border.color = node_colors[genotype[i]]
            face.margin_right = 0
            face.margin_left = 0
            face.border.width = 1

            face.inner_border.width = 0
            n.add_face(face, column=index, position="aligned")
            index += 1


    genotype_samples = {}
    leaves = set()
    for sample_id in leaf_meta:
        leaves.add(sample_id)
        if leaf_meta[sample_id][0] in genotype_samples:
            continue
        genotype_samples[leaf_meta[sample_id][0]] = sample_id

    valid_leaves = list(genotype_samples.values())
    ete_tree_obj.prune(valid_leaves)
    ete_tree_obj.render(outfile,tree_style=ts,dpi=dpi)
