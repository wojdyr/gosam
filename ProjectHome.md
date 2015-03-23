gosam is a bunch of Python scripts that can:
  * create monocrystal in PBC box (atomistic model),
  * create bicrystals with coincidence site lattice (CSL) boundaries in PBC box,
  * create crystalline grains of given shape, with vacancies, thermal vibrations, etc.
  * read/write several file formats (AtomEye cfg, VASP POSCAR, LAMMPS, DL\_POLY, XMOL XYZ).

**Dependencies:** Python and [numpy](http://numpy.scipy.org/). `graingen.py` also requires the [qhull](http://www.qhull.org/) program.

**Installation:** Just download the content of http://gosam.googlecode.com/svn/trunk/, preferably using SVN client.

**Documentation:** First read [Examples](Examples.md) and, if you are interested in generating grains, also [graingen](graingen.md). Most of the programs print basic usage instructions when invoked without parameters.

**See also:** [debyer](http://code.google.com/p/debyer/) - a bunch of tools coded in C and C++

contact: wojdyr@gmail.com