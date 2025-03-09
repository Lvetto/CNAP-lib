from ase.build import bulk

size = 10

# Gold (Au) Cubic Lattices
sc_au = bulk('Au', 'sc', a=4.08) * size     # Simple Cubic
bcc_au = bulk('Au', 'bcc', a=3.53) * size   # Body-Centered Cubic
fcc_au = bulk('Au', 'fcc', a=4.08) * size   # Face-Centered Cubic

# Save to XYZ files
sc_au.write(fr'data/sc{size}.xyz')
bcc_au.write(fr'data/bcc{size}.xyz')
fcc_au.write(fr'data/fcc{size}.xyz')

print("XYZ files for cubic gold lattices saved!")
