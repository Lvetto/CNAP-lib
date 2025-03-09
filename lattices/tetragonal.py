from ase import Atoms
from ase.io import write

# Define lattice parameters
a_tet = 4.08  # Å (assumed from cubic)
c_tet = 5.50  # Å (arbitrary)

size = 10

# Simple Tetragonal (1 atom per unit cell at (0,0,0))
st_au = Atoms(
    symbols="Au",
    scaled_positions=[[0, 0, 0]],  # Basis atom
    cell=[(a_tet, 0, 0), (0, a_tet, 0), (0, 0, c_tet)],  # Unit cell vectors
    pbc=True,  # Periodic boundary conditions
) * size

# Body-Centered Tetragonal (Extra atom at (a/2, a/2, c/2))
bct_au = Atoms(
    symbols="Au2",
    scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],  # Two atoms per unit cell
    cell=[(a_tet, 0, 0), (0, a_tet, 0), (0, 0, c_tet)],
    pbc=True,
) * size

# Save to XYZ file
write(fr'data/bct{size}.xyz', bct_au)
print("Body-Centered Tetragonal (BCT) saved as gold_bct.xyz")

# Save to XYZ file
write(fr'data/st{size}.xyz', st_au)
print("Simple Tetragonal (ST) saved as gold_st.xyz")
