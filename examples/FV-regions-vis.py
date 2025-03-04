from crystal_maker import *
import numpy as np
from CNA import *
from Graphs.FigureBuilder import *
from Graphs.Plots import *
from Graphs.DataClasses import *

# create a periodic lattice
a = 1
lattice = np.array(cubic_lattice_from_cell(1, 5, fcc_cell(a)))

# make a graph to represent the system
g = Graph(adjacency_matrix(0, lattice, 0.8*a), node_positions=lattice)

# extract subgraphs for different regions
fvs = get_feature_vectors(g)
keys, grouped_indices = group_fvs(fvs)
grouped_subgraphs = [g.make_subgraph(indices) for indices in grouped_indices]

# draw plots for each unique structure and its distribution in 3d space
fig2 = [FigureBuilder((1,1)) for key in keys]
for n, graph, sig in zip(range(len(grouped_subgraphs)), grouped_subgraphs, keys):
    fig2[n].add_plot(graph_plot, ax=fig2[n][0], graph=graph, pointsize=300, linesize=1, linealpha=0.6, pointalpha=0.8, title=sig)

plt.show()
