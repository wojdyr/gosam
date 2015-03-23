gosam is a bunch of Python scripts:
  * `bicryst.py` -- to generate atomistic model of grain boundary (bicrystal)
    with a coincidence site lattice (CSL) boundary in the PBC box,
  * `monocryst.py` -- to generate monocrystal in PBC box (atomistic model),
  * `graingen.py` -- create crystalline grains of given shape, with vacancies,
    thermal vibrations, etc.
  * `mdfile.py` -- convert between several file formats (AtomEye cfg,
    VASP POSCAR, LAMMPS, DL_POLY, XMOL XYZ).
  * `csl.py` -- to list coincident site lattice (CSL) grain boundaries
    in cubic crystals

**Prerequisites:** Python 2 and [numpy](http://numpy.scipy.org/).
`graingen.py` additionally requires the [qhull](http://www.qhull.org/) program.

**Usage**:
See [documentation in the wiki](https://github.com/wojdyr/gosam/wiki)
and feel free to e-mail me with questions.  Or open Github issue.
Also, most of the programs print basic usage instructions when invoked
without parameters.

**License**: GPLv2+ (I'll consider changing it to a more liberal license
if there is a good reason)

**Citing**:
you may cite
[Marcin Wojdyr et al 2010 Modelling Simul. Mater. Sci. Eng. 18 075009](
http://dx.doi.org/10.1088/0965-0393/18/7/075009)
([preprint](http://wojdyr.github.io/Wojdyr-tilt_GB_in_SiC-MSMSE-2010.pdf)),
especially if you use `bicryst.py`.

**See also:** [debyer](https://github.com/wojdyr/debyer) - a bunch of tools
coded in C and C++

contact: wojdyr@gmail.com
