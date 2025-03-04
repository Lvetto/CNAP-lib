from Graphs.FigureBuilder import *
from Graphs.Plots import *
from crystal_maker import *
from Graphs.DataClasses import *
import numpy as np
from CNA import *
from time import time
from math import ceil, pi

# create a periodic lattice
a = 1
lattice = np.array(cubic_lattice_from_cell(a, 15, fcc_cell(a)))

# get the signatures
adj_mat = adjacency_matrix(0, lattice, 0.8*a)
np.fill_diagonal(adj_mat, 0)
signatures = get_signature(adj_mat)

# create feature vectors
fvs = get_feature_vectors(signatures)

# find the geometric center of the system and compute distances from the center
center = np.sum(lattice, axis=0) / len(lattice)
radii = np.linalg.norm(lattice - center, axis=1)

# compile a list of unique local environments
unique_LAEs = []
for LAE in fvs.values():
    if (LAE not in unique_LAEs):
        unique_LAEs.append(LAE)

# make a set of bins and find the appropriate bin for each distance
r_min = np.min(radii)
r_max = np.max(radii)
dr = 0.5
bin_edges = np.array([r_min + i * dr for i in range(ceil((r_max - r_min) / dr) + 1)])
binned_radii = np.digitize(radii, bin_edges)
bins = [[0 for _ in bin_edges] for _ in unique_LAEs]

# make a list to keep track of occurrences of every LAE type
LAE_occurrences = [0 for _ in enumerate(unique_LAEs)]

# count LAEs for each bin
for atom, LAE in zip(fvs.keys(), fvs.values()):
    bins[unique_LAEs.index(LAE)][binned_radii[atom]] += 1
    LAE_occurrences[unique_LAEs.index(LAE)] +=1

volumes = [(4/3)*pi*((r+dr)**3-r**3) for r in bin_edges]    # volume of each spherical shell
densities = [[c/v for c,v in zip(bins[n], volumes)] for n,_ in enumerate(bins)]

# sort bar indices by their max height
sorted_indices = np.argsort(np.array([max(bin) for bin in bins]))

# make a figure and plot data on it
fig = FigureBuilder((1, 1), title="Radial distribution of LAEs")
for idx in reversed(sorted_indices):
    fig[0].bar(bin_edges, densities[idx], width=dr, label=f"{unique_LAEs[idx]}", edgecolor="k", alpha=0.5)

fig[0].legend()
fig[0].set_xlabel("distance from center")
fig[0].set_ylabel("density")

# make another figure to just show the number of occurrences of every LAE (I need the figure in the thesis)
fig2 = FigureBuilder((1,1), title="LAEs occurrences")
fig2[0].bar(range(len(LAE_occurrences)), LAE_occurrences, width=0.5, edgecolor="k")
fig2[0].set_xticks(range(len(LAE_occurrences)))
fig2[0].set_xticklabels(unique_LAEs)

plt.tight_layout()
plt.show()
