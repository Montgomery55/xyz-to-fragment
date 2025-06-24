import numpy as np
from scipy.spatial import distance_matrix
import networkx as nx
import pandas as pd
from itertools import permutations
from modules.bond_distances import *
from .vdw_radii import *
import time

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

    def ff_atom_types_basic(self, pd_df):
        ff_atom_types = []
        df = pd.read_csv(pd_df)
        for i, atom in enumerate(self.atoms):
            filtered = df[(df['# bonds']==self.bond_orders[i]) & (df['elemental symbol']==atom)]
            ff_atom_type = filtered['ff atom number'].iloc[0]
            ff_atom_types.append(ff_atom_type)
        return ff_atom_types

    def tinker_input_generator(self, df):
        N = np.arange(self.num_atoms)
        ff_atom_types = self.ff_atom_types_basic(df)
        print(self.num_atoms)
        for atom_index in N:
            entry = f'{atom_index+1:<4} {self.atoms[atom_index]:<8} {self.coords[atom_index][0]:>12.8f} {self.coords[atom_index][1]:>12.8f} {self.coords[atom_index][2]:>12.8f}\t{ff_atom_types[atom_index]:<6}\t{" ".join(map(str, 1+np.array(self.connectivities[atom_index])))}'
            print(entry)

    def ghost_atom_generator(self):
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
        fragmented_complex = np.array([np.array(x) for x in fragmented_complex], dtype=object)
        num_frags = np.arange(fragmented_complex.shape[0])

        print('-- dimer')
        for i in num_frags:
            for atom_coord in fragmented_complex[i]:
                print(f'{atom_coord[0]} {atom_coord[1]} {atom_coord[2]} {atom_coord[3]}')
        for i in num_frags:
            print(f'-- monomer {i+1}')
            for atom_coord in fragmented_complex[i]:
                print(f'{atom_coord[0]} {atom_coord[1]} {atom_coord[2]} {atom_coord[3]}')
            for j in num_frags:
                if j != i:
                    for atom_coord in fragmented_complex[j]:
                        print(f'@{atom_coord[0]} {atom_coord[1]} {atom_coord[2]} {atom_coord[3]}')
            
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

    def vdw_surface_area(self, grid_spacing=0.1):
        fragments = self.fragment()
        surface_areas = {}
        surface_points = {}
        for frag_num, frag in enumerate(fragments):
            coords = np.array([atom[1:] for atom in frag])
            atoms = np.array([atom[0] for atom in frag])
            radii = np.array([vdw_radii[atom[0]] for atom in frag])
            min_corner = coords.min(axis=0) - radii.max()
            max_corner = coords.max(axis=0) + radii.max()

            x = np.arange(min_corner[0], max_corner[0]+grid_spacing, grid_spacing)
            y = np.arange(min_corner[1], max_corner[1]+grid_spacing, grid_spacing)
            z = np.arange(min_corner[2], max_corner[2]+grid_spacing, grid_spacing)
            nx, ny, nz = len(x), len(y), len(z)

            occupancy = np.zeros((nx, ny, nz))
            occupancy = occupancy.astype(bool)
            xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')
            grid_points = np.stack([xv, yv, zv], axis=-1).reshape(-1, 3)

            grid_shape = (nx, ny, nz)
            grid_indices_3d = np.array(np.unravel_index(np.arange(grid_points.shape[0]), grid_shape)).T  # shape (n_points, 3)

            for coord, radius in zip(coords, radii):
                d2 = np.sum((grid_points - coord)**2, axis=1)
                inside = d2 <= radius**2
                occupancy_indices = grid_indices_3d[inside]

                occupancy[occupancy_indices[:, 0], occupancy_indices[:, 1], occupancy_indices[:, 2]] = True

            surface_area = 0.0
            face_area = grid_spacing**2
            directions = [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]
            surf_coords = []

            for dx, dy, dz in directions:
                shifted = np.roll(occupancy, shift=(dx, dy, dz), axis=(0, 1, 2))
                mask = occupancy & (~shifted)

                if dx == -1:
                    mask[0,:,:] = False
                elif dx == 1:
                    mask[-1,:,:] = False
                if dy == -1:
                    mask[:,0,:] = False
                elif dy == 1:
                    mask[:,-1,:] = False
                if dz == -1:
                    mask[:,:,0] = False
                elif dz == 1:
                    mask[:,:,-1] = False

                surface_area += np.sum(mask) * face_area
                surf_indices = np.argwhere(mask)
                surf_xyz = np.stack([x[surf_indices[:,0]],
                                        y[surf_indices[:,1]],
                                        z[surf_indices[:,2]]], axis=-1)
                surf_coords.append(surf_xyz)

            surface_coords = np.vstack(surf_coords)
            surface_areas[f'Frag {frag_num}'] = np.round(surface_area, 2)
            surface_points[f'Frag {frag_num}'] = surface_coords
        return surface_areas, surface_points







