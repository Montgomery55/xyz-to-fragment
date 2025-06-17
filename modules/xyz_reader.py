import numpy as np

class XYZ():
    def __init__(self, file):
        self.coords = None
        self.atoms = None
        self.file = file
        self.num_atoms = None

    def reader(self):
        num_atoms = 0
        atoms_coords = []
        with open(self.file, 'r') as f:
            for i, line in enumerate(f):
                if i == 0:
                    num_atoms = int(line)
                elif i > 1:
                    atoms_coords.append(line.split())
        atoms_coords = np.array(atoms_coords)
        atoms = atoms_coords[:,0]
        coords = atoms_coords[:,1:]
        coords = np.array([[float(val[0]), float(val[1]), float(val[2])] for val in coords])
        self.coords = coords
        self.atoms = atoms
        self.num_atoms = num_atoms

