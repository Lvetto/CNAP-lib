from ase.build import bulk
from ase.io import write

# Define lattice parameters for hexagonal gold
a_hex = 4.08  # Arbitrary choice
c_hex = 6.66  # c/a ratio ≈ 1.63 (ideal HCP value)

size = 10

# Create simple hexagonal structure
hex_au = bulk('Au', 'hcp', a=a_hex, c=c_hex) * size

# Save to XYZ file
write(fr'data/hex{size}.xyz', hex_au)
print("Simple Hexagonal (Hex) saved as gold_hex.xyz")
