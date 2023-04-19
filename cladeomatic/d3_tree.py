# Import library
from d3graph import d3graph, vec2adjmat

# Set source and target nodes
source = ['A','B','C','B','B','C']
target = ['A','A','A','D','E','F']
weight = [1,1,1,1,1,1]

# Create adjacency matrix
adjmat = vec2adjmat(source, target, weight=weight)

# target  node A  node B  node F  node J  node M  node C  node Z
# source
# node A    0.00     0.0    5.56    0.00    3.28     0.0     0.0
# node B    0.00     0.0    1.13    0.64    0.00     0.0     0.0
# node F    0.00     0.5    0.00    0.00    0.00     0.0     0.0
# node J    0.00     0.0    0.00    0.00    0.00     0.0     0.0
# node M    0.00     0.0    0.00    0.00    0.00     0.0     0.0
# node C    0.00     0.0    0.00    0.00    0.50     0.0     0.0
# node Z    0.45     0.0    0.00    0.00    0.00     0.0     0.0

# Initialize
d3 = d3graph()

# Build force-directed graph with default settings
d3.graph(adjmat)
d3.show(filepath='/Users/jrobertson/Desktop/temp.html')
