from Graphs.FigureBuilder import *
from Graphs.Plots import *
from crystal_maker import *
from Graphs.DataClasses import *
import numpy as np
from CNA import *
from math import ceil, sqrt
from random import random

# create a list of points
a = 1
cell = fcc_cell(a)

cutoff = 0.8 * a
adj_mat = adjacency_matrix(0, cell, cutoff)
bonds = get_all_bonds(adj_mat)

# make a list of random colors to make the plots easier to read
colors = [(random(), random(), random()) for _ in enumerate(cell)]

# prepare data to draw it
lines = [(cell[i], cell[j]) for i,j in bonds]
linesets = [LineSet(np.array(lines), alpha=1, linestyle=":")]
pointclouds = [PointCloud(np.array(cell), size = 1000, alpha=0.8, color=colors)]
labelsets = [LabelSet(np.array(cell), range(len(cell)))]

# draw the system
fig1 = FigureBuilder((1, 1))
fig1.add_plot(better_particle_plot, ax=fig1[0], points=pointclouds, lines=linesets, labels=labelsets, hide_frames=True)

# prepare the second figure with slots for each LAE
n = ceil(sqrt(len(cell)))
fig2 = FigureBuilder((n, n))

# get the adj list
adj_list = get_adj_list(adj_mat)

for atom, neighbors in enumerate(adj_list):
    atoms = [atom] + list(neighbors)
    positions = [cell[i] for i in atoms]

    lines = []
    for i in atoms:
        for j in adj_list[i]:
            if (j in atoms):
                lines.append((cell[i], cell[j]))

    # prepare data to plot
    pointclouds = [PointCloud(np.array(positions), size=300, alpha=0.8, color=[colors[i] for i in atoms])]
    labelsets = [LabelSet(np.array(positions), atoms)]
    linesets = [LineSet(np.array(lines), alpha=1, linestyle=":")]

    # plot data
    fig2.add_plot(better_particle_plot, ax=fig2[atom//n, atom%n], points=pointclouds, lines=linesets, labels=labelsets, hide_frames=True, title=f"{atom}")


plt.show()

