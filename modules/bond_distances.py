import numpy as np

covalent_radii_pm = {
    "H": 37,  "He": 28, #hydrogen gets a little fudge factor here, should be 31 but that seems too small
    "Li": 128, "Be": 96,  "B": 84,   "C": 77,   "N": 75,   "O": 66,   "F": 64,   "Ne": 58,
    "Na": 166, "Mg": 141, "Al": 121, "Si": 111, "P": 115,  "S": 105,  "Cl": 102, "Ar": 106,
    "K": 203,  "Ca": 176, "Sc": 170, "Ti": 160, "V": 153,  "Cr": 139, "Mn": 139, "Fe": 132,
    "Co": 126, "Ni": 124, "Cu": 132, "Zn": 122, "Ga": 122, "Ge": 120, "As": 119, "Se": 120,
    "Br": 120, "Kr": 116,
    "Rb": 220, "Sr": 195, "Y": 190,  "Zr": 175, "Nb": 164, "Mo": 154, "Tc": 147, "Ru": 146,
    "Rh": 142, "Pd": 139, "Ag": 145, "Cd": 144, "In": 142, "Sn": 139, "Sb": 139, "Te": 138,
    "I": 139,  "Xe": 140,
    "Cs": 244, "Ba": 215, "La": 207, "Ce": 204, "Pr": 203, "Nd": 201, "Pm": 199, "Sm": 198,
    "Eu": 198, "Gd": 196, "Tb": 194, "Dy": 192, "Ho": 192, "Er": 189, "Tm": 190, "Yb": 187,
    "Lu": 187,
    "Hf": 175, "Ta": 170, "W": 162,  "Re": 151, "Os": 144, "Ir": 141, "Pt": 136, "Au": 136,
    "Hg": 132, "Tl": 145, "Pb": 146, "Bi": 148, "Po": 140, "At": 150, "Rn": 150,
    "Fr": 260, "Ra": 221, "Ac": 215, "Th": 206, "Pa": 200, "U": 196, "Np": 190, "Pu": 187,
    "Am": 180, "Cm": 169
}

def bonding_radius(atom, charge=0, ion_factor=0.70):
    if charge < 0: #anions are larger than neutrals
        return covalent_radii_pm[atom]/ion_factor * 0.01 #0.01 converts pm to angstroms
    elif charge > 0: #cations are smaller than neutrals
        return covalent_radii_pm[atom]*ion_factor * 0.01 #0.01 converts pm to angstroms
    else:
        return covalent_radii_pm[atom] * 0.01 #0.01 converts pm to angstroms


