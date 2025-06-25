# Fragmentation Library

This is a library for obtaining noncovalent fragments of xyz files capable of producing tinker input files, ghost atom fragmentation, and general fragment xyz files.

## Usage

Input sould be an `.xyz` file with the following format:

```
# of atoms

atom1 x1 y1 z1
atom2 x2 y2 z2
...
atomN xN yN zN
```

The `modules.xyz_reader.XYZ` class provides several main methods:
`tinker_input_generator('file')`: produces a tinker input file from a given `.xyz` file.

`fragment()`: produces a fragment (noncovalent) of a given `.xyz` file.

`ghost_atom_generator`: produces an `.xyz` file where ghost atoms are given as '@X'.

`vdw_surface_area(grid_spacing=n)`: produces both the van der Waal surface area and the surface coordinates. This returns a tuple of `surface_area` in Angstrom$^2$ and `surface_points`. 
The density of the points is controled by `grid_spacing` (a value of 0.2 is sufficient in most cases, but the accuracy can be increased by decreasing the value of `grid_spaxing`, but if the system is large, 0.5 will work).
The size of the surface is controlled by `solvent_probe` which adds a value to the van der Waal radii and/pr `vdw_scale` which scales the van der Waal radii by a scalar value (between 0 and 1).

`ESP_map()`: produces normalized Gasteiger charges (obtained from the electronegativity and chemical hardness (calculated via IP and EA values) at a given surface (such as that obtained from `vdw_surfave_area`).

These core methods are dependent on two other methods of the `XYZ` class, `reader()` and `bond_order_connectivities()`, which must be run before the three other methods given.

The `main.py` provides examples on how to run these methods.

## Installation

Run `pip install -r requirements.txt` in your python environment to install dependencies.

## Tests

To confirm that the code is running correctly, a `test_main.py` file is given.
Running the test file, `python test_main.py`, will compare the given outputs to what they should be.
If the test(s) are past, one will see `file.out passes` printed to the terminal.

