from crystal_maker import *
import numpy as np
from CNA import *
from Graphs.FigureBuilder import *
from Graphs.Plots import *
from Graphs.DataClasses import *
from random import random
from old.pp_io import read_xyz

# read data from file
filepath = r"data/bcc5.xyz"
elements, points = read_xyz(filepath)

cutoff = 1.1

# make a graph to represent the system
g = Graph(adjacency_matrix(0, points, cutoff), node_positions=points)

# create feature vectors, group them and format the data
fvs = get_feature_vectors(g)
keys, grouped = group_fvs(fvs)

print(f"Unique LAEs found: {keys}")

# some config for the figures
fps=8
rot_steps = (0, 4, 0)
rot0 = (20, 0, 0)

# draw a figure showing the spatial distribution of FVs
fig = FigureBuilder((1, 1))
colors = [(random(), random(), random()) for _ in grouped]
pointclouds = [PointCloud(points[i], size=200, color=j, label=k, alpha=1) for i,j,k in zip(grouped, colors, keys)]
fig.add_plot(rotating_particle_plot, plot_pos=0, points=pointclouds, hide_frames=False, title="", legend=True, rot_steps=rot_steps, rot0=rot0)
fig.show(frames=90, save=True, path=rf"data/out/all.gif", fps=fps)
del fig
print(f"saved figure: all")

# save figures for different FVs
for key, inds, color in zip(keys, grouped, colors):
    fig = FigureBuilder((1, 2), title=key)
    fig.add_plot(rotating_graph_plot, plot_pos=0, graph=make_neighbors_subgraph(g, inds[0]), rot_steps=rot_steps, rot0=rot0 ,hide_frames=True, pointsize=100, linesize=1, linealpha=0.6, pointalpha=0.8)
    fig.add_plot(rotating_graph_plot, plot_pos=1, graph=g.make_subgraph(inds), rot_steps=rot_steps, rot0=rot0 ,hide_frames=True, pointsize=50, linesize=1, linealpha=0.6, pointalpha=0.8, pointcolor=color)
    
    fig.show(frames=90, save=True, path=rf"data/out/{key}.gif", fps=fps)
    del fig
    print(f"saved figure: {key}")
