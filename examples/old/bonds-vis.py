from Graphs.FigureBuilder import *
from Graphs.Plots import *
from crystal_maker import *
from Graphs.DataClasses import *
import numpy as np
from CNA import *

# create the figure object
fig = FigureBuilder((1, 3))

# create primitive cells for the 3 cubic lattices
a = 1
cells = [sc_cell(a), bcc_cell(a), fcc_cell(a)]
titles = ["Simple cubic", "Body centered cubic", "Face centered cubic"]
cutoffs = [1.1 * a, 1.2 * a, 0.8 * a]

# compute adj matrix and extract all bonds
adj = [adjacency_matrix(0, i[0], i[1]) for i in zip(cells, cutoffs)]
bonds = [get_all_bonds(i) for i in adj]

# create data needed to draw lines to represent the bonds
# probably could be made into a list comprehension, but doesn't sound like an enjoyable endevour
lines = []
for cell_bonds, cell_points in zip(bonds, cells):
    cell_lines = [(cell_points[i], cell_points[j]) for i,j in cell_bonds]
    lines.append(np.array(cell_lines))

# prepare data for the plot
linesets = [LineSet(i, alpha=1, linestyle=":") for i in lines]
pointclouds = [PointCloud(np.array(i), size = 1000, alpha=0.8) for i in cells]

# draw the cubic cells
[fig.add_plot(better_particle_plot, ax=fig[i], points=[pointclouds[i]], lines=[j[1]], hide_frames=True, title=j[0]) for i,j in enumerate(zip(titles, linesets))]

plt.show()
