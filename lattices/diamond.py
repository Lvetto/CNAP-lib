from ase.lattice.cubic import Diamond
from ase.io import write
from ase.neighborlist import NeighborList
import numpy as np

# Define the element and lattice constant
element = "Au"  # Gold
a = 4.08  # Lattice constant for diamond structure in Å

size = 5

# Generate diamond lattice
atoms = Diamond(symbol=element, latticeconstant=a, size=(size, size, size))

# Save as XYZ file
write("data/diamond10.xyz", atoms)

# Compute first-neighbor distance using NeighborList
cutoff = 0.6 * a  # Set a reasonable cutoff for neighbors
nl = NeighborList([cutoff] * len(atoms), self_interaction=False, bothways=True)
nl.update(atoms)

# Get first neighbor distances for the first atom
indices, offsets = nl.get_neighbors(0)
distances = [atoms.get_distance(0, i, mic=True) for i in indices]

# First neighbor distance is the smallest
first_neighbor_distance = min(distances)

# Print result
print(f"First neighbor distance for the diamond lattice: {first_neighbor_distance:.4f} Å")
