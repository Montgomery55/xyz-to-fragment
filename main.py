from modules.xyz_reader import *
from modules.bond_distances import *

file = 'tests/water_cluster.xyz'

sample = XYZ(file)
sample.reader()
sample.bond_order_connectivities()

sample.tinker_input_generator('modules/basic_ff_atom_types.txt')
print()
sample_frag = sample.fragment()
for frag in sample_frag:
    for atom in frag:
        print(atom[0], atom[1], atom[2], atom[3])
    print()
print()
sample.ghost_atom_generator()
