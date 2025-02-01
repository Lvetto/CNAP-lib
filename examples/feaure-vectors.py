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

# get all unique signatures
unique_sigs = list(set(signatures.values()))

# create feature vectors
fv = get_feature_vectors(signatures)

# find the geometric center of the system and compute distances from the center
center = np.sum(lattice, axis=0) / len(lattice)
radii = np.linalg.norm(lattice - center, axis=1)

# make a set of bins and find the appropriate bin for each distance
r_min = np.min(radii)
r_max = np.max(radii)
dr = 0.3
bin_edges = np.array([r_min + i * dr for i in range(ceil((r_max - r_min) / dr) + 1)])
binned_radii = np.digitize(radii, bin_edges)
bins = [[0 for _ in bin_edges] for _ in unique_sigs]

# count signatures in each bin
for atom, sigs in zip(fv.keys(), fv.values()):
    for sig, count in zip(sigs.keys(), sigs.values()):
        bins[unique_sigs.index(sig)][binned_radii[atom]] += count

# make a figure and fill it with bar plots
fig = FigureBuilder((1, len(unique_sigs) + 1), title="Number of atoms participating in a signature")
for n, sig in enumerate(unique_sigs):
    t = [bin_edges, bins[n]]
    fig.add_plot(bar_plot, ax=fig[n], bar_data=t, title=f"{unique_sigs[n]}", width=dr, x_title="Radius", y_title="Number of atoms")

# compute the volume of each spherical shell and plot them
volumes = [(4/3)*pi*((r+dr)**3-r**3) for r in bin_edges]
fig[-1].plot(volumes)
fig[-1].set_title("V(r+dr)-V(r)")

# make another figure and plot the densities
fig2 = FigureBuilder((1, len(unique_sigs)), title="Density of atoms participating in a signature")
for n, sig in enumerate(unique_sigs):
    t = [bins[n][i]/volumes[i] for i,_ in enumerate(volumes)]
    t = [bin_edges, t]
    fig2.add_plot(bar_plot, ax=fig2[n], bar_data=t, title=f"{unique_sigs[n]}", width=dr, x_title="Radius", y_title="Density of atoms")

plt.show()
