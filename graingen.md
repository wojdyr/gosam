# graingen.py #

### Installation ###

#### Debian/Ubuntu ####
```
$ sudo apt-get subversion python-numpy qhull-bin
$ svn checkout http://gosam.googlecode.com/svn/trunk/ gosam
```

#### Windows ####

Install Python and Numpy. Use subversion client (e.g. TortoiseSVN) to download _gosam_. Download the [qhull](http://www.qhull.org/) program and copy the file qhull.exe to the `gosam` directory.

#### Other systems ####
Install Python, Numpy and qhull. Use subversion client to download _gosam_.

### Basic usage ###

`graingen.py` generates crystalline grains. It is recommended to use it from a python script that looks like a configuration file, but has imports at the beginning and a call to `generate_grain()` at the end.

Example 1. Let's create a file Fe\_r4.py:

```
from graingen import *

cell = CubicUnitCell(a=2.87)

atoms = [
    ("Fe", 0.0, 0.0, 0.0),
    ("Fe", 0.5, 0.5, 0.5),
]

surfaces = LatticeSurface(r=4)

generate_grain(globals())
```

Put this file in the `gosam` directory.

On Windows: click the script to execute it.

On Unix system: type
```
$ python Fe_r4.py
```

The script should produce two files with the same basename: log file `Fe_r4.log` and output configuration `Fe_r4.xyz`.

![http://gosam.googlecode.com/svn/wiki/Fe-r4.png](http://gosam.googlecode.com/svn/wiki/Fe-r4.png)

The import statement above, in the first line of the script, works
if your script is in the same directory as the `graingen.py` file.
Alternatively, you may install _gosam_ as Python package and use:
```
from gosam.graingen import *
```

Assignment of the `cell` variable defines unit cell. Usually values in ångströms are used, but the program is unit agnostic. Any unit can be used (consistently). The coordinates in the output file will be in the same unit that was used in the configuration file.

Assignment of the `atoms` variable defines positions of atoms in the unit cell.

The `surfaces` variable is responsible for size and shape of the grain. In the example above, the surface is spherical and the radius is 4Å.

The script above generates bcc lattice and cuts out an atomic cluster discarding atoms outside of the r=4Å sphere.

Instead of the sphere, one can use a convex polyhedron defined by a set of planes:

```
  surfaces = [
    # each plane is defined by Miller indices and distance r from (0,0,0)
    LatticeSurface( hkl=(0,-1,0), r=25 ), 
    LatticeSurface( hkl=(-1,0,0), r=25 ),
    # ...
    # polyhedron needs at least 4 planes
  ]
```

### Surface deformation ###

It is possible to define surface deformation (the relaxation of atoms at the surface). The deformation can be defined as displacement of atoms from the perfect lattice positions, in direction normal to the surface, with displacement magnitude given as a function of the distance from the surface.
For example:
```
LatticeSurface.default_sd = SurfaceDeformation(depth=8.3, 
                                               fun=lambda x: (8.3-x)*0.05 )
```
To assign different deformation for each surface plane define LatticeSurface with parameter `sd`, like this:
```
def func(d):
    return 0.02 * (10-d)**2

surfaces = [
    LatticeSurface( hkl=(1,0,0),  r=5 ),
    LatticeSurface( hkl=(-1,0,0), r=5 ),
    LatticeSurface( hkl=(0,1,0),  r=10 ),
    LatticeSurface( hkl=(0,-1,0), r=10 ),
    LatticeSurface( hkl=(0,0,1),  r=15, sd=SurfaceDeformation(10, func) ),
    LatticeSurface( hkl=(0,0,-1), r=15 ),
]
```

![http://gosam.googlecode.com/svn/wiki/deform.png](http://gosam.googlecode.com/svn/wiki/deform.png)

### Thermal vibrations ###

If the function named `modifier` is defined, it is used to modify each atom.
This can be used to add thermal vibrations:
```
def modifier(atom):
    sigma = 0.025
    atom.pos += (random.gauss(0, sigma),
                 random.gauss(0, sigma),
                 random.gauss(0, sigma))
```

Each atom has three properties:
  * `name` -- atomic symbol
  * `pos` -- position (x, y, z) in ångströms
  * `min_dist` -- distance to the grain surface
It is possible to specify different vibrations for atoms near the surface, or different vibrations for different atomic species.

### Vacancies ###

It is possible to generate vacancies (remove some atoms) by specifying the probability of vacancy:
```
vacancy_probability = 0.01
```

Different species can have different probabilities:
```
vacancy_probability = { "Si": 0.01, "C": 0.02 }
```

The `vacancy_probability` variable can be defined as a function that takes atom and returns probability:
```
def vacancy_probability(atom):
    if atom.min_dist < 5: # surface, 5A shell
        return 0.05
    else                  # core
        return 0.01
```

### Atom substitutions ###

The modifier() function can be used to add substitutional impurity atoms:
```
# this substitutes, with 2% probability, any atoms with oxygen
def modifier(atom):
    if random.random() < 0.02:
        atom.name = "O"
```

Note that if you already have the modifier() function in your file, you need to add code to that function. Do not create two functions with the same name.

### Nodes ###

Instead of `atoms`, as in the first example, two variables: `nodes` and `node_atoms` can be used.

```
nodes = [
    (0.0, 0.0, 0.0),
    (0.5, 0.5, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
]

node_atoms = [
    ("Si", 0.0, 0.0, 0.0),
    ("C",  0.25,0.25,0.25),
]
```

There are two predefined lists of nodes: `fcc_nodes` and `bcc_nodes`. This is equivalent to the assignment in the example above:
```
nodes = fcc_nodes
```

If the unit cell has _n_ nodes and _m_ node atoms, it has _n_ × _m_ atoms.
The positions in the second list are added to the node position.

Nodes contain one or more atoms that are, by default, kept together. To ungroup nodes, add line:
```
do_not_group_nodes=True
```

### Output ###

The `output_file` variable sets the basename of log file and output files.
If not set, the basename is set to the basename of the script.
```
output_file = "foo"
```

The `output_formats` variable contains a list of output formats. It should contain one or more formats supported by `mdfile.py` (XYZ, AtomEye, LAMMPS and others).
```
output_formats = ["xyz"]
output_formats = ["xyz", "atomeye"]
```

Additionally, a log of the operation is written to a file with the `.log` extension.

### Example ###

More complex example that generates a spherical grain of cubic SiC:

```
#!/usr/bin/env python

import random
from graingen import *

output_file = "cub_sic"
output_formats = ["atomeye"]

cell = CubicUnitCell(4.36)

#nodes in unit cell (as fraction of unit cell parameters)
nodes = [
    (0.0, 0.0, 0.0),
    (0.5, 0.5, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
]
# alternatively, predefined lists fcc_nodes and bcc_nodes can be used, e.g:
#nodes = fcc_nodes

#atoms in node (as fraction of unit cell parameters)
node_atoms = [
    ("Si", 0.0, 0.0, 0.0),
    ("C",  0.25,0.25,0.25),
]

# the simplest case - spherical grain
surfaces = LatticeSurface(r=25)

# add thermal vibrations
def modifier(atom):
    sigma = 0.025
    atom.pos += (random.gauss(0, sigma),
                 random.gauss(0, sigma),
                 random.gauss(0, sigma))

#vacancy_probability = { "Si": 0.01, "C": 0.02 }

generate_grain(globals())
```

![http://gosam.googlecode.com/svn/wiki/cub_sic.png](http://gosam.googlecode.com/svn/wiki/cub_sic.png)