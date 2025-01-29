from Graphs.FigureBuilder import *
from Graphs.Plots import *
from crystal_maker import *
from Graphs.DataClasses import *
import numpy as np
from CNA import *
from time import time

t0 = time()

a = 1
lattice = np.array(cubic_lattice_from_cell(a, 15, fcc_cell(a)))

t1 = time()

adj_mat = adjacency_matrix(0, lattice, 0.8*a)
np.fill_diagonal(adj_mat, 0)
signatures = get_signature(adj_mat)

t2 = time()

occurances = {}
for bond, signature in zip(signatures.keys(), signatures.values()):
    if (signature in occurances.keys()):
        occurances[signature] += 1
    else:
        occurances[signature] = 1

t3 = time()

print(f"Number of atoms {len(lattice)}\nTime taken:\n\tCreating lattice: {t1-t0:.2f}\n\tComputing signatures: {t2-t1:.2f}\n\tCounting occurances: {t3-t2:.2f}")
print("Occurances:")
[print(f"\t{signature}: {occurrance}") for signature, occurrance in zip(occurances.keys(), occurances.values())]
