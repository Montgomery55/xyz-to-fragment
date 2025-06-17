import numpy as np
from scipy.spatial import distance_matrix
import networkx as nx
import pandas as pd
from itertools import permutations
from modules.bond_distances import *

class XYZ():
    def __init__(self, file):
        self.coords = None
        self.atoms = None
        self.file = file
        self.num_atoms = None
        self.bond_orders = None
        self.connectivities = None

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

    def bond_matrix(self):
        dist_matrix = distance_matrix(self.coords, self.coords)
        N = self.coords.shape[0]
        bonded_matrix = np.zeros(dist_matrix.shape)
        for i, j in np.ndindex(bonded_matrix.shape):
            distance = dist_matrix[i][j]
            bond_i = bonding_radius(self.atoms[i])
            bond_j = bonding_radius(self.atoms[j])
            bond_max = bond_i + bond_j
            if distance <= bond_max and distance != 0:
                bonded_matrix[i][j] = 1
        return bonded_matrix

    def bond_order_connectivities(self):
        bond_order = []
        bonding_matrix = self.bond_matrix()
        for atom in bonding_matrix:
            bond_order.append(int(sum(atom)))
        connectivities = [list(np.where(row==1)[0]) for row in bonding_matrix]
        self.bond_orders = np.array(bond_order)
        self.connectivities = connectivities

    def ff_atom_types(self, pd_df):
        ff_atom_types = []
        df = pd.read_csv(pd_df)
        for i, atom in enumerate(self.atoms):
            filtered = df[(df['# bonds']==self.bond_orders[i]) & (df['elemental symbol']==atom)]
            ff_atom_type = filtered['ff atom number'].iloc[0]
            ff_atom_types.append(ff_atom_type)
        return ff_atom_types

    def tinker_input_generator(self, df):
        N = np.arange(self.num_atoms)
        self.ff_atom_types(df)
        ff_atom_types = self.ff_atom_types(df)
        print(self.num_atoms)
        for atom_index in N:
            entry = f'{atom_index+1:<4} {self.atoms[atom_index]:<8} {self.coords[atom_index][0]:>12.8f} {self.coords[atom_index][1]:>12.8f} {self.coords[atom_index][2]:>12.8f}\t{ff_atom_types[atom_index]:<6}\t{" ".join(map(str, 1+np.array(self.connectivities[atom_index])))}'
            print(entry)


    def fragment(self):
        bonding_matrix = self.bond_matrix()
        G = nx.Graph()
        G.add_nodes_from(range(self.num_atoms))

        for i in range(self.num_atoms):
            for j in range(self.num_atoms):
                if bonding_matrix[i, j] == 1:
                    G.add_edge(i, j)

        fragments = list(nx.connected_components(G))
        fragmented_complex = []
        for frags in fragments:
            fragment_n = []
            for i in frags:
                fragment_n.append([self.atoms[i], self.coords[i][0], self.coords[i][1], self.coords[i][2]])
            fragmented_complex.append(fragment_n)
        return fragmented_complex

