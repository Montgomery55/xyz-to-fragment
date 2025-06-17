from modules.xyz_reader import *
from modules.bond_distances import *
import io
from contextlib import redirect_stdout

test_complicated_molecule = 'tests/complicated_molecule.xyz'
test_benzene = 'tests/benzene.xyz'
test_methane_dimer = 'tests/methane_dimer.xyz'
test_water_cluster = 'tests/water_cluster.xyz'

complicated_molecule = XYZ(test_complicated_molecule)
complicated_molecule.reader()
complicated_molecule.bond_order_connectivities()
f = io.StringIO()
with redirect_stdout(f):
    complicated_molecule.tinker_input_generator('modules/basic_ff_atom_types.txt')
output = f.getvalue()
with open('tests/complicated_molecule.out', 'r') as f:
    expected = f.read()
if output == expected:
    print('complicated_molecule.xyz passes')
else:
    print('something is wrong with output file for complicated_molecule.xyz')

benzene = XYZ(test_benzene)
benzene.reader()
benzene.bond_order_connectivities()
f = io.StringIO()
with redirect_stdout(f):
    benzene.tinker_input_generator('modules/basic_ff_atom_types.txt')
output = f.getvalue()
with open('tests/benzene.out', 'r') as f:
    expected = f.read()
if output == expected:
    print('benzene.xyz passes')
else:
    print('something is wrong with output file for benzene.xyz')

methane_dimer = XYZ(test_methane_dimer)
methane_dimer.reader()
methane_dimer.bond_order_connectivities()
f = io.StringIO()
with redirect_stdout(f):
    methane_dimer.tinker_input_generator('modules/basic_ff_atom_types.txt')
output = f.getvalue()
with open('tests/methane_dimer.out', 'r') as f:
    expected = f.read()
if output == expected:
    print('methane_dimer.xyz passes')
else:
    print('something is wrong with output file for methane_dimer.xyz')

water_cluster = XYZ(test_water_cluster)
water_cluster.reader()
water_cluster.bond_order_connectivities()
f = io.StringIO()
with redirect_stdout(f):
    water_cluster.tinker_input_generator('modules/basic_ff_atom_types.txt')
output = f.getvalue()
with open('tests/water_cluster.out', 'r') as f:
    expected = f.read()
if output == expected:
    print('water_cluster.xyz passes')
else:
    print('something is wrong with output file for water_cluster.xyz')
