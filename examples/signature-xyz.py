from Graphs.FigureBuilder import *
from Graphs.Plots import *
from crystal_maker import *
from Graphs.DataClasses import *
import numpy as np
from CNA import *
from time import time
from old.pp_io import read_xyz
from scipy.spatial import ConvexHull

t0 = time()

# read data from file
filepath = r"data/tail.xyz"
elements, points = read_xyz(filepath)

# make a convex hull approximation for the volume of the system and compute the number density
hull = ConvexHull(points)
volume = hull.volume
density = len(points) / volume

# make a rough extimate for the cutoff distance as proportional to the average distance between points
alpha = 1.2
cutoff = alpha * (density)**(-1/3)

t1 = time()

# get the signatures
adj_mat = adjacency_matrix(0, points, cutoff)
np.fill_diagonal(adj_mat, 0)
signatures = get_signature(adj_mat)

t2 = time()

# count occurrances
occurances = {}
for bond, signature in zip(signatures.keys(), signatures.values()):
    if (signature in occurances.keys()):
        occurances[signature] += 1
    else:
        occurances[signature] = 1

t3 = time()

# print out some info
print(f"Number of atoms {len(points)}\nTime taken:\n\tCreating lattice: {t1-t0:.2f}s\n\tComputing signatures: {t2-t1:.2f}s\n\tCounting occurances: {t3-t2:.2f}s")
print("Occurances:")
[print(f"\t{signature}: {occurrance}") for signature, occurrance in zip(occurances.keys(), occurances.values())]

# draw the points and their bonds
pointsets = [PointCloud(points, size=10, alpha=0.5)]
bonds = get_all_bonds(adj_mat)
line_data = np.array([(points[i], points[j]) for i,j in bonds])
linesets = [LineSet(line_data, linestyle=":", alpha=1)]

fig = FigureBuilder((1, 1))
fig.add_plot(better_particle_plot, ax=fig[0], points=pointsets, lines=linesets, hide_frames=False)
plt.show()
