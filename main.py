from modules.xyz_reader import *
from modules.bond_distances import *

file = 'tests/water_cluster.xyz'

benzene = XYZ(file)
benzene.reader()
"""
fragmented_complex = benzene.fragment()
for frag in fragmented_complex:
    for atom in frag:
        print(atom[0], atom[1], atom[2], atom[3])
    print()
"""
benzene.bond_order_connectivities()

coords = benzene.coords
atoms = benzene.atoms
atom_num = benzene.num_atoms
bond_orders = benzene.bond_orders
atom_connectivities = benzene.connectivities

print(bond_orders)
print(atom_connectivities)
print(benzene.bond_matrix())

