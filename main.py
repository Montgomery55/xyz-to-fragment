from modules.xyz_reader import *
from modules.bond_distances import *

file = 'tests/benzene.xyz'

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

benzene.tinker_input_generator('modules/basic_ff_atom_types.txt')
