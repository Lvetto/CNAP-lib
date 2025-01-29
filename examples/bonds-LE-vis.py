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
np.fill_diagonal(adj_mat, 0)
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
n = ceil(sqrt(len(bonds)))
fig2 = FigureBuilder((n, n))

adj_list = get_adj_list(adj_mat)
CNs = get_common_neighbors(adj_list, bonds)

for count, pair, neighbors in zip(range(len(CNs.keys())), CNs.keys(), CNs.values()):
    atoms = list(pair) + list(neighbors)
    positions = [cell[i] for i in atoms]

    bonds = get_common_bonds(adj_mat, neighbors)

    lines = []
    for i,j in bonds:
        lines.append((cell[i], cell[j]))

    # prepare data to plot
    pointclouds = [PointCloud(np.array(positions), size=300, alpha=0.8, color=[colors[i] for i in atoms])]
    labelsets = [LabelSet(np.array(positions), atoms)]
    linesets = [LineSet(np.array(lines), alpha=1, linestyle=":")]

    # put the two atoms considered and the first two values of the CNA signature as a title
    title = f"Bond {pair[0]}-{pair[1]} ({len(neighbors)} {len(bonds)})"

    # plot data
    fig2.add_plot(better_particle_plot, ax=fig2[count//n, count%n], points=pointclouds, lines=linesets, labels=labelsets, hide_frames=True, title=title)


plt.show()
