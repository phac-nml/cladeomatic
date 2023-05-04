import json

import matplotlib
import networkx as nx
import pandas as pd
import seaborn as sns
from deprecated import deprecated

@deprecated
def create_clade_heirarchy(genotypes,delim="."):

    h = {}
    for sample_id in genotypes:
        genotype = genotypes[sample_id].split(delim)
        if len(genotype) <= 1:
            continue
        for i in range(1,len(genotype)):
            h[genotype[i]] = genotype[i-1]

    return h
@deprecated
def create_graph_json(h,outfile):
    g = nx.DiGraph()

    for child in h:
        parent = h[child]
        g.add_edge(parent, child)
    root = [n for n, d in g.in_degree() if d == 0][0]

    cmap = sns.color_palette('husl', n_colors=len(g.nodes()))
    count = 0
    for node in nx.dfs_preorder_nodes(g, root):
        g.nodes[node]['color'] = matplotlib.colors.rgb2hex(cmap[count])
        count += 1


    tree = nx.tree_data(g, root, ident="name")
    json.dump(tree, open(outfile, 'w'), indent=2)


genotypes = {
    'root':'0',
    'root-1':'0.1',
    'root-2':'0.2',
    'root-1.1': '0.1.1',
    'root-1.2': '0.1.2',
    'root-2.1': '0.2.1',
    'root-2.2': '0.2.1',
}