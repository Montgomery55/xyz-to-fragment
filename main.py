from modules.xyz_reader import *

file = 'tests/benzene.xyz'

benzene = XYZ(file)
benzene.reader()
coords = benzene.coords
atoms = benzene.atoms
atom_num = benzene.num_atoms

print(coords)
print(atoms)
print(atom_num)
