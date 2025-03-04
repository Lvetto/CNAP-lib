from crystal_maker import *
import numpy as np
from CNA import *
from Graphs.FigureBuilder import *
from Graphs.Plots import *
from Graphs.DataClasses import *
from random import random

# create a periodic lattice
a = 1
lattice = np.array(cubic_lattice_from_cell(1, 5, bcc_cell(a)))

# make a graph to represent the system
g = Graph(adjacency_matrix(0, lattice, 1.1*a), node_positions=lattice)

# create feature vectors, group them and format the data
fvs = get_feature_vectors(g)
keys, grouped = group_fvs(fvs)

# create figures
fig = FigureBuilder((1, 1))
fig2 = [FigureBuilder((1,2), title=key) for key in keys]#FigureBuilder((2, len(keys)), title="LAEs")

# some data for the plots
colors = [(random(), random(), random()) for _ in grouped]
pointclouds = [PointCloud(lattice[i], size=300, color=j, label=k, alpha=1) for i,j,k in zip(grouped, colors, keys)]

# add subplots and show them
fig.add_plot(better_particle_plot, plot_pos=0, points=pointclouds, hide_frames=False, title="", legend=True)

# draw plots for each unique structure and its distribution in 3d space
for n, key, inds in zip(range(len(grouped)), keys, grouped):
    fig2[n].add_plot(graph_plot, plot_pos=0, graph=make_neighbors_subgraph(g, inds[0]), pointsize=300, linesize=1, linealpha=0.6, pointalpha=1, title=key)
    fig2[n].add_plot(graph_plot, plot_pos=1, graph=g.make_subgraph(inds), pointsize=300, linesize=1, linealpha=0.6, pointalpha=1, title=key, pointcolor=colors[n])

plt.show()
