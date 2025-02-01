from Graphs.FigureBuilder import *
from Graphs.Plots import *
from crystal_maker import *
from Graphs.DataClasses import *
import numpy as np
from CNA import *
from time import time

t0 = time()

# create a periodic lattice
a = 1
lattice = np.array(cubic_lattice_from_cell(a, 10, fcc_cell(a)))

t1 = time()

# get the signatures
adj_mat = adjacency_matrix(0, lattice, 0.8*a)
np.fill_diagonal(adj_mat, 0)
signatures = get_signature(adj_mat)

t2 = time()

# count occurrances
occurrences = get_occurrences(signatures)

t3 = time()

# Print out some info
print(f"Number of atoms {len(lattice)}\nTime taken:\n\tCreating lattice: {t1-t0:.2f}s\n\tComputing signatures: {t2-t1:.2f}s\n\tCounting occurrences: {t3-t2:.2f}s")
print("Occurrences:")
[print(f"\t{signature}: {occurrance}") for signature, occurrance in zip(occurrences.keys(), occurrences.values())]
