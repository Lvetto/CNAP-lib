from crystal_maker import *
import numpy as np
from CNA import *
from Graphs.FigureBuilder import *
from Graphs.Plots import *
from Graphs.DataClasses import *
from random import random

# create a periodic lattice
a = 1
lattice = np.array(cubic_lattice_from_cell(1, 3, fcc_cell(a)))

# make a graph to represent the system
g = Graph(adjacency_matrix(0, lattice, 0.8*a), node_positions=lattice)

# create feature vectors, group them and format the data
fvs = get_feature_vectors(g)
keys, grouped = group_fvs(fvs)
positions = [np.array([i.position for i in group]) for group in grouped]

# create figures
fig = FigureBuilder((1, 2))
fig2 = [FigureBuilder((1,2), title=key) for key in keys]#FigureBuilder((2, len(keys)), title="LAEs")

# some data for the plots
colors = [(random(), random(), random()) for _ in positions]
pointclouds = [PointCloud(i, size=100, color=j, label=k, alpha=0.8) for i,j,k in zip(positions, colors, keys)]

# add subplots and show them
fig.add_plot(graph_plot, ax=fig[0], graph=g, pointsize=100, linesize=1, linealpha=0.6, pointalpha=0.8)
fig.add_plot(better_particle_plot, ax=fig[1], points=pointclouds, hide_frames=False, title="", legend=True)

# draw plots for each unique structure and its distribution in 3d space
for n, key in enumerate(keys):
    i = find_index(fvs, "signatures", key)
    sub_g = make_neighbors_subgraph(g, i)
    fig2[n].add_plot(graph_plot, ax=fig2[n][0], graph=sub_g, pointsize=300, linesize=1, linealpha=0.6, pointalpha=0.8, title=key)
    fig2[n].add_plot(better_particle_plot, ax=fig2[n][1], points=[pointclouds[n]], hide_frames=False, title="", legend=True)


fig.show()
