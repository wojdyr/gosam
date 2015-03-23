# mdfile.py #
Converts between file formats. Filenames with `gz` and `bz2` extensions are recognized and handled as expected. See the list of supported files [in the source](http://code.google.com/p/gosam/source/browse/trunk/mdfile.py#4).

Convert VASP POSCAR file to LAMMPS input.
```
mdfile.py POSCAR foo.lammps./monocryst.py  --margin=1 NaCl 4 4 4 nacl.cfg
```

Convert all AtomEye cfg files to gzipped input LAMMPS files.
```
bash$ for i in *.cfg; do ../../mdfile.py $i `basename $i .cfg`.lammps.gz; done
```

Change all atoms B to C and Cl to Si.
```
mdfile.py --translate 'B->C, Cl->Si' input.cfg output.cfg
```

Swap Si and C atoms.
```
mdfile.py --translate 'Si->C, C->Si' input.cfg output.cfg
```

Set explicitly PBC box size.
```
mdfile.py --pbc '[(60,0,0),(0,60,0),(0,0,60)]' input.xyz output.cfg
```

Remove atoms that have the z coordinate outside of the (15,30) range.
```
mdfile.py --filter '15 < z < 30' input.cfg output.cfg
```


---


# csl.py #
Info about coincident site lattice (CSL) grain boundaries in cubic crystals.

Show usage instructions.
```
csl.py
```

Print angle of rotation for Σ=5 grain boundary with the `[`001] axis of rotation.
```
csl.py 001 5
```


---


# pse.py #
Periodic System of Elements

```
$ python -c 'import pse; print pse.pse_dict["Au"]'
Z=79 gold (Au)  mass=196.967
```


---


# monocryst.py #
Generates monocrystal in PBC. Out of box supports Fe (bcc), Cu (fcc), NaCl, Si and diamond (both have the diamond structure), SiC (zinc blende, wurtzite or any polytype specified by a sequence of letters). Other systems can be easily added in the [get\_named\_lattice()](http://code.google.com/p/gosam/source/browse/trunk/monocryst.py#188) function.

Generate Si crystal with the size equal or larger than 2nm x 2nm x 3nm.
```
monocryst.py si 2 2 3 output.cfg
```

![http://gosam.googlecode.com/svn/wiki/Si-223.png](http://gosam.googlecode.com/svn/wiki/Si-223.png)

Option `--margin` generates a cuboid in PBC, by padding the monocrystal with vacuum (of course it's not a monocrystal anymore).
```
monocryst.py --margin=1 NaCl 4 4 4 nacl.cfg
```
![http://gosam.googlecode.com/svn/wiki/nacl-cubes.png](http://gosam.googlecode.com/svn/wiki/nacl-cubes.png)

All the three pictures present the same configuration, shifted under PBC.
The command above generates the system in the left picture.
With the additional `--center-zero` option, the center of the configuration is exactly at (0,0,0), like in the middle picture.
In AtomEye program, it is possible to shift the view under PBC, using mouse. If your visualization program cannot do this and you want to see the cube in one piece, you may use a program from the [debyer](http://code.google.com/p/debyer/) toolset:
```
# shift the system by half of the PBC in x, y and z directions.
dbr_extend -r -S0.5,0.5,0.5 -i nacl.cfg
```

Having a configuration with the center at or near (0,0,0) makes mathematic calculations easier. For example, command
```
mdfile.py --prefer-negative --filter='abs(x)+abs(y)+abs(z) < 20' nacl.cfg nacl-oct.cfg
```
makes octahedron. The `--prefer-negative` option changes the internal processing of the coordinates: they are mapped to the (-H/2, H/2) range, instead of (0, H), where H is the PBC box size.

![http://gosam.googlecode.com/svn/wiki/nacl-oct.png](http://gosam.googlecode.com/svn/wiki/nacl-oct.png)



---


# bicrystal.py #
Generates bicrystals.
This script has been used to generate only symmetric tilt and twist GBs in cubic systems.
Nonetheless, it has a lot of options.

Print help (but the help is not complete).
```
bicrystal.py
```

Generate diamond bicrystal:
  * rotation axis `[`100],
  * median plane (011), i.e. its a tilt boundary,
  * Σ=13 (pick the lowest rotation angle for Σ=13),
  * min. dimensions: 0.8 x 0.8 x 3 nm,
  * do not adjust z dimension
  * add 1.0nm vacuum in z dimension
  * if two atoms are in a distance smaller than 0.5, remove one of them
```
bicrystal.py 100 m011 13 0.8 0.8 3 nozfit remove:0.5 vacuum:1.0 lattice:diamond d13-ini.cfg
```

The GB is constructed at the z=0 plane, so you may want to shift it under PBC before visualization, like in the nacl.cfg case above.

![http://gosam.googlecode.com/svn/wiki/d13-ini.png](http://gosam.googlecode.com/svn/wiki/d13-ini.png)



---


Let's now try to construct `<`110> symmetric tilt grain boundaries in Cu, similar to that in the paper: M. Tschopp _et al_., Acta Materialia 55 (2007) 3959-3969.
Let's pick Σ267 (11,11,5) θ=144.4° boundary (Fig. 3f).
To learn more about `<`110> Σ267 GBs, use csl.py:

```
$ ./csl.py 110 267
m=22  n= 5  35.636
m=13  n= 7  74.579
m=14  n=13 105.421
m= 5  n=11 144.364
```

The m and n integers are used to generate boundaries (see Grimmer, Acta Cryst. (1984) A40, 108). m,n can be used as a bicrystal.py parameter, instead of the Σ.
In this case: `5,11`.

If the rotation axis is `[`110], the (11,11,5) plane does not contain the axis, so to have tilt boundary, we use the equivalent (11,-11,5) plane.

The final command is
```
./bicrystal.py 1,1,0 11,-11,5 5,11 2 2 3 lattice:cu out.cfg
```
And we got GB with the 5 macroscopic degrees of freedom the same as in the Fig. 3f.

![http://gosam.googlecode.com/svn/wiki/cu267.png](http://gosam.googlecode.com/svn/wiki/cu267.png)

The system is periodic and there are two GBs. To get only one GB use the `vacuum` option, as in the previous example.

As you can see, in this geometrical construction some GB atoms are very close to each other. This system is clearly not a realistic model of the GB. The `remove` option (see the previous example) can remove the atoms that are too close to each other, but this is only the first step (and this step may not be necessary) to create realistic models.

The complete procedure that we used to prepare GBs is described in the paper:<br>
M. Wojdyr <i>et al</i>, <i>Energetics and structure of <code>&lt;001&gt;</code> tilt grain boundaries in SiC</i>,<br>
<a href='http://dx.doi.org/10.1088/0965-0393/18/7/075009'>Modelling Simul. Mater. Sci. Eng. 18 075009 (2010)</a>.<br>
<br>
<br>
<hr />

<h1>graingen.py</h1>

<a href='graingen.md'>described here</a>

<hr />
<h1>extra features</h1>

Branden Kappes has added a method <code>orthogonalize_pbc()</code> to convert models from nonorthogonal to orthogonal PBCs. It does not lead to periodic boundaries but rather creates a fully dense, orthogonal box. His use case was:<br>
<br>
<blockquote>I'm working with graphene (two triangular lattices offset by (a+b)/3 where a, b are the periodicity vectors--graphene, being 2D, does not have a meaningful c vector) on various substrates.</blockquote>

<blockquote>(1) Because of its triangular structure, the most logical place to start is graphene on (111) surfaces of metals (Ni, Ir, etc.) that closely match the periodicity of the graphene. (The periodicity of the Ni(111) triangular lattice is only ~1.5% larger than that of graphene.)</blockquote>

<blockquote>(2) Because I am using DFT as implemented in VASP, I need to minimize the size of my systems. The most expedient way to do that is by using non-orthogonal models that match the underlying periodicity of the constituent crystal systems (alpha, beta, gamma = 60, 120, 90 for graphene on FCC(111)).  VASP handles these just fine, but the ReaxFF work I want to do in LAMMPS does not.