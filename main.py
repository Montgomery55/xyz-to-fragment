from modules.xyz_reader import *
from modules.bond_distances import *
import numpy as np

file = 'tests/water_cluster.xyz'

sample = XYZ(file)
sample.reader()
sample.bond_order_connectivities()

#sample tinker input generator
sample.tinker_input_generator('modules/basic_ff_atom_types.txt')
print()

#sample fragmentation
sample_frag = sample.fragment()
for frag in sample_frag:
    for atom in frag:
        print(atom[0], atom[1], atom[2], atom[3])
    print()
print()

#sample ghost atom generator
sample.ghost_atom_generator()
print()

#sample surface area and surface coordinate generator
surface_area, surface_points = sample.vdw_surface_area(grid_spacing=0.2)
print(f'van der Waal surface area of water (fragment 1): {np.round(surface_area["Frag 0"], 2)} \u212b\u00b2')
print(f'Surface coordinates of Fragment 1 {surface_points["Frag 0"]}')

