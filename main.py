from modules.xyz_reader import *
from modules.bond_distances import *

file = 'tests/MeNH2_peptide.xyz'

test = XYZ(file)
test.reader()
test.bond_order_connectivities()

test.ghost_atom_generator()
