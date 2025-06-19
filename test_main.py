from modules.xyz_reader import *
from modules.bond_distances import *
import io
from contextlib import redirect_stdout

test_complicated_molecule = 'tests/complicated_molecule.xyz'
test_benzene = 'tests/benzene.xyz'
test_methane_dimer = 'tests/methane_dimer.xyz'
test_water_cluster = 'tests/water_cluster.xyz'
test_MeNH2_peptide = 'tests/MeNH2_peptide.xyz'

print('----------------------------------')
print('FF INPUT FILE GENERATION TESTS')
print('----------------------------------')
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
    print('complicated_molecule.out passes')
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
    print('benzene.out passes')
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
    print('methane_dimer.out passes')
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
    print('water_cluster.out passes')
else:
    print('something is wrong with output file for water_cluster.xyz')

print('----------------------------------')
print('FRAGMENTATION TEST')
print('----------------------------------')
fragmented_complex = water_cluster.fragment()
f = io.StringIO()
with redirect_stdout(f):
    for frag in fragmented_complex:
        for atom in frag:
            print(atom[0], atom[1], atom[2], atom[3])
        print()
output = f.getvalue()
with open('tests/water_cluster_frag.out', 'r') as f:
    expected = f.read()
if output == expected:
    print('water_cluster_frag.out passes')
else:
    print('something is wrong with output file for water_cluster_frag.xyz')

print('----------------------------------')
print('GHOST FUNCTION TESTS')
print('----------------------------------')
MeNH2_peptide = XYZ(test_MeNH2_peptide)
MeNH2_peptide.reader()
MeNH2_peptide.bond_order_connectivities()
f = io.StringIO()
with redirect_stdout(f):
    MeNH2_peptide.ghost_atom_generator()
output = f.getvalue()
with open('tests/MeNH2_peptide_ghosts.out', 'r') as f:
    expected = f.read()
if output == expected:
    print('MeNH2_peptide.out passes')
else:
    print('something is wrong with output file for MeNH2_peptide.xyz')

water_cluster = XYZ(test_water_cluster)
water_cluster.reader()
water_cluster.bond_order_connectivities()
f = io.StringIO()
with redirect_stdout(f):
    water_cluster.ghost_atom_generator()
output = f.getvalue()
with open('tests/water_cluster_ghosts.out', 'r') as f:
    expected = f.read()
if output == expected:
    print('water_cluster.out passes')
else:
    print('something is wrong with output file for water_cluster.xyz')
