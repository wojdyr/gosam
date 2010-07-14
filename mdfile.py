#!/usr/bin/env python
# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2
"""
File import/export.
Supported file types (some of them are supported partially):
 XMOL XYZ File Format (search WWW for details)
 DL_POLY (see http://www.cse.clrc.ac.uk/msi/software/DL_POLY/)
 AtomEye (search WWW for details)
 LAMMPS (http://lammps.sandia.gov)
 Pielaszek (simple file format used sometimes here, in Unipress)
 POSCAR - VASP input file with atom positions - POSCAR (direct format)
 GULP - a part of GULP input file with cell vectors and atom coordinates
"""

import sys
import operator
from optparse import OptionParser
import bz2
import gzip
import numpy
from numpy import linalg

from mdprim import AtomVF
import model
import pse
from utils import get_command_line
from rotmat import is_diagonal, StdDev


def get_orthorhombic_pbc(full_pbc):
    """\
    the input is matrix 3x3, the output is diagonal.
    Non-diagonal elements must be zero.
    """
    if full_pbc is None or len(full_pbc) == 0:
        raise ValueError("no PBC")
    assert is_diagonal(numpy.array(full_pbc))
    return numpy.diagonal(full_pbc)


def export_for_dlpoly(atoms, f, title, sort=True):
    levcfg = 0 #only coordinates in file
    #imcon = 0 #no periodic boundaries
    imcon = 1 #cubic periodic boundaries
    a = 200.
    with_shells = False
    print >>f, title
    print >>f, "%10i%10i" % (levcfg, imcon)
    print >>f, "%20f%20f%20f" % (a, 0, 0)
    print >>f, "%20f%20f%20f" % (0, a, 0)
    print >>f, "%20f%20f%20f" % (0, 0, a)
    ordered = atoms[:]
    if sort:
        # atoms are sorted in alphabetical order
        ordered.sort(key=operator.attrgetter("name"))

    # output atoms to file
    counter = 0
    for atm in ordered:
        counter += 1
        print >>f, "%-8s%20i" % (atm.name, counter)
        print >>f, "%20f%20f%20f" % (atm.pos[0], atm.pos[1], atm.pos[2])
        if with_shells: #shell - hack
            counter += 1
            print >>f, "%-8s%20i" % (atm.name + "-shl", counter)
            print >>f, "%20f%20f%20f" % (atm.pos[0], atm.pos[1], atm.pos[2])

    # gather data useful for FIELD file
    atom_counts = {}
    for i in ordered:
        if i.name not in atom_counts:
            atom_counts[i.name] = 1
        else:
            atom_counts[i.name] += 1
    return atom_counts


def export_as_xmol(atoms, f, title):
    print >>f, len(atoms)
    print >>f, title
    for point in atoms:
        print >>f, point.name, point.pos[0], point.pos[1], point.pos[2]


def export_for_pielaszek(atoms, f):
    "the simplest format, used sometimes here in Unipress"
    for point in atoms:
        print >>f, point.pos[0], point.pos[1], point.pos[2], point.name


def import_pielaszek(ifile):
    "the simplest format, used sometimes here in Unipress"
    atoms = []
    for i in ifile:
        s = i.split()
        pos = (float(s[0]), float(s[1]), float(s[2]))
        atoms.append(AtomVF(s[3], len(atoms), pos, None, None))
    return model.Model(atoms, pbc=[])


def import_xmol(ifile):
    number_of_atoms = int(ifile.readline())
    title = ifile.readline().strip()
    atoms = []
    for i in ifile:
        s = i.split()
        pos = (float(s[1]), float(s[2]), float(s[3]))
        atoms.append(AtomVF(s[0], len(atoms), pos, None, None))
    return model.Model(atoms, pbc=[], title=title)


def import_dlpoly_config(ifile):
    title = ifile.readline().strip()
    second_line = ifile.readline().split()
    levcfg = int(second_line[0])
    imcon = int(second_line[1])
    return _get_dlpoly_configuration(ifile, title, levcfg, imcon)[0]


def _get_dlpoly_configuration(ifile, title, levcfg, imcon):
    def get_numbers(line):
        return [float(i) for i in line.split()]
    assert imcon <= 3, "Sorry, no support for non-parallelepiped PBC"
    pbc = []
    if imcon > 0:
        for i in range(3):
            pbc.append(get_numbers(ifile.readline()))
    atoms = []
    next_line = None
    while 1:
        record1 = ifile.readline()
        if not record1:
            break
        if record1.startswith("timestep"):
            next_line = record1
            break
        rec1_sp = record1.split()
        name = rec1_sp[0]
        nr = int(rec1_sp[1])
        pos = get_numbers(ifile.readline())
        velocity = None
        force = None
        if levcfg > 0:
            velocity = get_numbers(ifile.readline())
            if levcfg > 1:
                force = get_numbers(ifile.readline())
        atoms.append(AtomVF(name, nr, pos, velocity, force))
    return model.Model(atoms, pbc, title), next_line


def import_dlpoly_history(ifile):
    title = ifile.readline().strip()
    second_line = ifile.readline().split()
    levcfg = int(second_line[0])
    imcon = int(second_line[1])
    all_configurations = []
    next_line = ifile.readline()
    while next_line:
        sys.stdout.write(".")
        sys.stdout.flush()
        ts = next_line.split()
        assert int(ts[3]) == levcfg and int(ts[4]) == imcon
        c, next_line = _get_dlpoly_configuration(ifile, title, levcfg, imcon)
        #TODO yield c
        all_configurations.append(c)
    sys.stdout.write("\n")
    return all_configurations


def dlpoly_history_info(ifile):
    title = ifile.readline().strip()
    print "title:", title
    second_line = ifile.readline().split()
    print "number of atoms in frame:", second_line[2]
    print "has following frames:",
    frame_counter = 0
    while 1:
        line = ifile.readline()
        if not line:
            break
        if line.startswith("timestep"):
            print line.split()[1],
            sys.stdout.flush()
            frame_counter += 1
    print "finished.", frame_counter, "frames were found."


def get_stoichiometry_string(configuration):
    counts = configuration.count_species()
    return "Stoichiometry: " + " ".join("%s:%d" % i for i in counts.iteritems())


def export_for_atomeye(configuration, f):
    "AtomEye Extended CFG format"
    pbc = configuration.pbc
    atoms = configuration.atoms
    aux = [
            # uncomment to add temperature (works if there are atom velocities
            # in the input file)

           #("temperature [K]", lambda a: a.get_temperature()),

            # add aux value that can be used to color code atomic positions
            # in unit cell, in selected direction. The coloring works with hsv
            # color scale in AtomEye, or any other scale that is cyclic, i.e.
            # the same color corresponds to values 0 and 1.

           #("xcolor [0-1]", in_cell_pos_fun(0, pbc[0][0]/2,
           #                                   pos0=_find_pos0(atoms))),
          ]

    # atoms from VASP POSCAR with Selective dynamics have allow_change attr
    if hasattr(atoms[0], "allow_change"):
        aux += [ ("change x", lambda a: float(a.allow_change[0])),
                 ("change y", lambda a: float(a.allow_change[1])),
                 ("change z", lambda a: float(a.allow_change[2])) ]

    if pbc is None or len(pbc) == 0:
        raise ValueError("no PBC")
    if not isinstance(pbc, numpy.ndarray):
        pbc = numpy.array(pbc)
    print >>f, "Number of particles = %i" % len(atoms)
    for i in get_comment_list(configuration):
        print >>f, "# " + i
    print >>f, "A = 1.0 Angstrom (basic length-scale)"
    for i in range(3):
        for j in range(3):
            print >>f, "H0(%i,%i) = %f A" % (i+1, j+1, configuration.pbc[i][j])
    print >>f, ".NO_VELOCITY."
    print >>f, "entry_count = %i" % (3 + len(aux))
    for n, a in enumerate(aux):
        print >>f, "auxiliary[%i] = %s" % (n, a[0])
    H_1 = linalg.inv(pbc)
    previous_name = None
    for i in atoms:
        if previous_name != i.name:
            print >>f, pse.get_atom_mass(i.name)
            print >>f, i.name
            previous_name = i.name
        s = numpy.dot(i.pos, H_1) % 1.0
        entries = [s[0], s[1], s[2]]
        for aname, afunc in aux:
            entries.append(afunc(i))
        print >>f, " ".join("%f" % i for i in entries)


# This function returns reference 0 point for coloring based on
# the x coordinate. The point is to color similar systems in the same way,
# even if the center of mass has been shifted.
# You may need to tune it for your data.
def _find_pos0(atoms):
    # in my systems B atoms are in a rigid surface, so we use it as
    # a reference
    selected = [a for a in atoms if a.name == "B"]
    if selected:
        elem = min(selected, key = lambda x: x[3])
        return elem[1]
    else: # no rigid surface
        #return sum(a.pos[0] for a in atoms) / len(atoms)
        return atoms[0].pos[0]


# returns function that returns value in the [0,1) range that corresponds to 
# position of atom in unit cell
def in_cell_pos_fun(dir, cell_size, pos0=0):
    print "adding aux for %s position in cell size: %g" % (chr(ord('x')+dir),
                                                           cell_size)
    return lambda a: ((a.pos[dir] - pos0) / cell_size) % 1.


def import_atomeye(ifile):
    # read header
    pbc = [[0,0,0], [0,0,0], [0,0,0]]
    has_velocity = True
    comments = ""
    for line in ifile:
        if not line:
            continue
        elif line[0] == '#':
            comments += line[1:]
        elif line.startswith("Number of particles"):
            number_of_atoms = int(line.split("=")[1])
        elif line.startswith("H0"):
            k,v = line[2:].split("=")
            a,b = k.strip().lstrip("(").rstrip(")").split(",")
            pbc[int(a)-1][int(b)-1] = float(v.split()[0])
        elif line.strip() == ".NO_VELOCITY.":
            has_velocity = False
        elif line[0].isdigit():
            break

    atoms = []
    vel = None
    H = numpy.array(pbc, float)
    for line in ifile:
        if not line or line.isspace() or line[0] == '#':
            continue
        elif line[0].isalpha():
            spec = line.strip()
        else:
            s = line.split()
            if len(s) == 1: # this is probably atomic mass
                continue
            pos = (float(s[0]), float(s[1]), float(s[2]))
            pos = numpy.dot(pos, H)
            if has_velocity:
                vel = (float(s[3]), float(s[4]), float(s[5]))
                # if velocities are also reduced:
                #vel = numpy.dot(vel, H)
            atoms.append(AtomVF(spec, len(atoms), pos, vel, None))
    return model.Model(atoms, pbc=pbc, title="from cfg", comments=comments)


# This file format doesn't contain atom names, only numbers of atom types.
# The names corresponding to the numbers are in separate lammps file.
# Here we assume that names of the types follow the line "n atom types"
# as comments, e.g.: "2 atom types # C Si"
def import_lammps_data(ifile):
    pbc = [[0,0,0], [0,0,0], [0,0,0]]
    comments = ""
    while 1:
        line = ifile.readline()
        if '#' in line:
            hash_idx = line.find('#')
            comment = line[hash_idx+1:]
            line = line[:hash_idx].strip()
        else:
            comment = ""
            line = line.strip()
        if not line:
            if comment:
                comments += comment
            continue
        if line.endswith("atoms"):
            number_of_atoms = int(line.split()[0])
        elif line.endswith("atom types"):
            number_of_atom_types = int(line.split()[0])
            atom_names = comment.split()
            assert len(atom_names) == number_of_atom_types
        elif line.endswith("xlo xhi"):
            lo, hi = line.split()[0:2]
            length = float(hi) - float(lo)
            pbc[0][0] = length
        elif line.endswith("ylo yhi"):
            lo, hi = line.split()[0:2]
            length = float(hi) - float(lo)
            pbc[1][1] = length
        elif line.endswith("zlo zhi"):
            lo, hi = line.split()[0:2]
            length = float(hi) - float(lo)
            pbc[2][2] = length
        elif line == "Atoms":
            break
        else:
            assert 0

    atoms = []
    while len(atoms) < number_of_atoms:
        s = ifile.readline().split()
        if len(s) < 5:
            continue
        n = int(s[0])
        assert n == len(atoms) + 1
        spec = int(s[1])
        name = atom_names[spec-1]
        pos = float(s[2]), float(s[3]), float(s[4])
        atoms.append(AtomVF(name, len(atoms), pos, vel=None, force=None))

    return model.Model(atoms, pbc=pbc, title="from lammps data",
                       comments=comments)

def get_comment_list(configuration):
    cmd = get_command_line()
    cc = ["file written by gosam (SVN $Revision$)", cmd]
    if configuration.title and configuration.title != cmd:
        cc.append(configuration.title)
    cc.append(get_stoichiometry_string(configuration))
    return cc


def export_as_lammps(configuration, f):
    "exporting coordinates as LAMMPS input data file"
    #configuration.orthogonalize_pbc()
    ort_pbc = get_orthorhombic_pbc(configuration.pbc)
    for i in get_comment_list(configuration):
        print >>f, "# " + i
    print >>f
    counts = configuration.count_species()
    species = sorted(counts.keys())
    print >>f, "\n%d\tatoms" % len(configuration.atoms)

    ## temporary hack: for now i need 4 atom types for 2 species
    #species += ["B", "Ge"]

    print >>f, "%d atom types # %s" % (len(species), " ".join(species))
    spmap = dict((i, n+1) for n, i in enumerate(species))
    print >>f, "0 %.6f xlo xhi" % ort_pbc[0]
    print >>f, "0 %.6f ylo yhi" % ort_pbc[1]
    print >>f, "0 %.6f zlo zhi" % ort_pbc[2]
    print >>f, "\nAtoms\n"
    for n, i in enumerate(configuration.atoms):
        print >>f, "%d\t%d\t%.7f\t%.7f\t%.7f" % (n+1, spmap[i.name],
                                           i.pos[0], i.pos[1], i.pos[2])


def export_as_poscar(configuration, f):
    "exporting coordinates as VASP POSCAR file"
    # make list of species
    counts = configuration.count_species()
    # reverse if you want to have Si before C
    species = sorted(counts.keys(), reverse=False)

    # line 1: a comment (usually the name of the system; we put here elements)
    print >>f, " ".join(species)
    # line 2: scaling factor
    print >>f, "%.15f" % 1
    # line 3, 4, 5: the unit cell of the system
    for i in range(3):
        print >>f, "%19.15f %19.15f %19.15f" % tuple(configuration.pbc[i])

    # line 6: the number of atoms per atomic species
    for i in species:
        print >>f, counts[i],
    print >>f

    # optional line
    selective_dynamics = hasattr(configuration.atoms[0], "allow_change")
    if selective_dynamics:
        print >>f, "Selective dynamics"

    # line 7: Cartesian/Direct switch
    print >>f, "Direct"

    # line 8, ...: atom coordinates
    H_1 = linalg.inv(configuration.pbc)
    for sp in species:
        for i in configuration.atoms:
            if i.name == sp:
                s = numpy.dot(i.pos, H_1) % 1.0
                if selective_dynamics:
                    allow_change = [('T' if i else 'F') for i in s[3:6]]
                    print >>f, "%19.15f %19.15f %19.15f %s %s %s" % (
                                s[0], s[1], s[2], c[0], c[1], c[2])
                else:
                    print >>f, "%19.15f %19.15f %19.15f" % tuple(s)


def import_poscar(ifile):
    first_line = ifile.readline().strip()
    # by our convention, first line contains symbols of elements
    species = first_line.split()
    scaling_factor = float(ifile.readline())
    assert scaling_factor > 0, \
           "POSCAR import: negative scaling factors are not supported"

    pbc = []
    for i in range(3):
        line = ifile.readline()
        h = [scaling_factor * float(i) for i in line.split()]
        assert len(h) == 3
        pbc.append(h)

    line = ifile.readline()
    atom_count = [int(i) for i in line.split()]

    # check if it's possible that the species above are set correctly
    assert len(species) == len(atom_count), \
            "POSCAR import: first line should contain (only) symbols of atoms"

    H = numpy.array(pbc, float)

    selective_dynamics = False
    cartesian = False
    switch = ifile.readline()[0].lower()
    if switch == 's':
        selective_dynamics = True
        switch = ifile.readline()[0].lower()
    if switch in ('c', 'k'):
        cartesian = True

    atoms = []
    for n, count in enumerate(atom_count):
        name = species[n]
        for i in xrange(count):
            s = ifile.readline().split()
            raw_pos = (float(s[0]), float(s[1]), float(s[2]))
            if cartesian:
                pos = scaling_factor*numpy.array(raw_pos)
            else:
                pos = numpy.dot(raw_pos, H)
            atom = AtomVF(name, len(atoms), pos, None, None);
            if len(s) == 6: # interpret T/F
                atom.allow_change = (s[3] != 'F', s[4] != 'F', s[5] != 'F')
            atoms.append(atom)

    return model.Model(atoms, pbc=pbc)


def export_as_gulp(configuration, f):
    "export coordinates in GULP input format"
    print >>f, ""
    print >>f, "# title:", configuration.title
    print >>f, "# Configuration"
    print >>f, "cell"
    pbc = get_orthorhombic_pbc(configuration.pbc)
    print >>f, "%.6g %.6g %.6g %g %g %g" % (pbc[0], pbc[1], pbc[2],
                                            90., 90., 90.)
    print >>f, "fractional"
    H_1 = linalg.inv(configuration.pbc)
    for i in configuration.atoms:
        s = numpy.dot(i.pos, H_1) % 1.0
        print >>f, "%s   core %11.6g %11.6g %11.6g    0.0000000  1.0000000" % (
                                                      i.name, s[0], s[1], s[2])

def get_type_from_filename(name):
    lname = name.lower()
    if name.endswith(".cfg"):
        return "atomeye"
    if name.endswith(".lammps") or name.endswith(".lmps"):
        return "lammps"
    elif name.endswith(".xyz"):
        return "xmol"
    elif name.endswith(".at"):
        return "pielaszek"
    elif name.endswith(".gin"):
        return "gulp"
    #TODO strip directory(?)
    elif "revcon" in lname or "config" in lname:
        return "dlpoly"
    elif "history" in lname:
        return "dlpoly_history"
    elif "poscar" in lname or "contcar" in lname:
        return "poscar"
    elif name.endswith(".bz2"):
        return get_type_from_filename(name[:-4])
    elif name.endswith(".gz"):
        return get_type_from_filename(name[:-3])
    else:
        print "Can't deduce filetype from filename:", name
        return None


def open_any(name, mode='r'):
    if name.endswith(".bz2"):
        return bz2.BZ2File(name, mode)
    elif name.endswith(".gz"):
        return gzip.GzipFile(name, mode)
    elif name == '-':
        if 'w' in mode:
            return sys.stdout
        else:
            return sys.stdin
    else:
        return file(name, mode)


def import_autodetected(filename):
    input_type = get_type_from_filename(filename)
    infile = open_any(filename)
    if input_type == "xmol":
        return import_xmol(infile)
    elif input_type == "dlpoly":
        return import_dlpoly_config(infile)
    elif input_type == "dlpoly_history":
        return import_dlpoly_history(infile)
    elif input_type == "pielaszek":
        return import_pielaszek(infile)
    elif input_type == "atomeye":
        return import_atomeye(infile)
    elif input_type == "lammps":
        return import_lammps_data(infile)
    elif input_type == "poscar":
        return import_poscar(infile)
    else:
        assert 0

def parse_translate_option(str):
    tr_map = {}
    for i in str.split(","):
        if "->" not in i:
            raise ValueError("wrong format of --translate option")
        k, v = i.split("->")
        tr_map[k.strip()] = v.strip()
    return tr_map


def parse_options(argv):
    parser = OptionParser("usage: %prog [--pbc=list] input_file output_file")
    parser.add_option("--pbc",
                      help="set PBC, e.g. '[(60,0,0),(0,60,0),(0,0,60)]'")
    parser.add_option("--center-zero", action="store_true",
                      help="move center of mass (with all atoms equal) to 0")
    parser.add_option("--prefer-negative", action="store_true",
                      help="if in PBC, select image in (-h/2, h/2)")
    parser.add_option("--filter",
                      help="criterion for leaving atoms, e.g. '15 < z < 30'")
    parser.add_option("--translate", metavar="TR",
                      help="change atomic symbols, e.g. 'Si->C, C->Si'")
    parser.add_option("--reference", help="adds dx, dy and dz properties",
                      metavar="FILE")
    (options, args) = parser.parse_args(argv)
    if len(args) == 0:
        parser.print_help()
        sys.exit()
    if not (len(args) == 2 or ("vs" in args and len(args) in (5,6))):
        parser.error("Two arguments (input and output filenames) are required")
    return (options, args)


def put_pbc_image_between_halfs(cfg):
    pbc = numpy.diagonal(cfg.pbc)
    for atom in cfg.atoms:
        for i in range(3):
            if atom.pos[i] > pbc[i] / 2.:
                atom.pos[i] -= pbc[i]
            elif atom.pos[i] < -pbc[i] / 2.:
                atom.pos[i] += pbc[i]


def process_input(input_filename, options):
    cfg = import_autodetected(input_filename)
    if options.pbc:
        cfg.pbc = eval(options.pbc)

    if options.prefer_negative:
        if not cfg.pbc:
            print "Error: Option --prefer-negative can be used only in PBC"
            sys.exit()
        put_pbc_image_between_halfs(cfg)

    if options.reference:
        cref = import_autodetected(options.reference)
        assert len(cref.atoms) == len(cfg.atoms)
        pbc = numpy.diagonal(cfg.pbc)
        for n, a1 in enumerate(cfg.atoms):
            a0 = cref.atoms[n]
            a1.dpos = [a1.pos[0] - a0.pos[0],
                       a1.pos[1] - a0.pos[1],
                       a1.pos[2] - a0.pos[2]]
            for i in range(3):
                if a1.dpos[i] > pbc[i] / 2.:
                    a1.dpos[i] -= pbc[i]
                elif a1.dpos[i] < -pbc[i] / 2.:
                    a1.dpos[i] += pbc[i]
        del cref

    if options.center_zero:
        ctr = sum(atom.pos for atom in cfg.atoms) / len(cfg.atoms)
        print "center: (%.3f, %.3f, %.3f)" % tuple(ctr)
        for atom in cfg.atoms:
            atom.pos -= ctr

    if options.filter:
        def f(atom):
            name = atom.name
            x, y, z = atom.pos
            return eval(options.filter)
        cfg.atoms = [i for i in cfg.atoms if f(i)]
        print len(cfg.atoms), "atoms left."
        if len(cfg.atoms) == 0:
            sys.exit()

    if options.translate:
        tr_map = parse_translate_option(options.translate)
        for i in cfg.atoms:
            if i.name in tr_map:
                i.name = tr_map[i.name]

    return cfg

def export_autodetected(configuration, output_filename):
    output_type = get_type_from_filename(output_filename)
    #TODO: check if file already exists
    ofile = open(output_filename, 'w')
    configuration.export_atoms(output_filename, output_type)

def convert(argv):
    "converts atomistic files"
    options, args = parse_options(argv)
    configuration = process_input(args[0], options)
    export_autodetected(configuration, args[1])


def get_atom_func(name):
    def x(atom): return atom.pos[0]
    def y(atom): return atom.pos[1]
    def z(atom): return atom.pos[2]
    def vx(atom): return atom.vel[0] # [A/ns]
    def vy(atom): return atom.vel[1]
    def vz(atom): return atom.vel[2]
    def v(atom): return atom.get_velocity()
    def T(atom): return atom.get_temperature()
    def Ekin(atom): return atom.get_ekin()
    def dx(atom): return atom.dpos[0]
    def dy(atom): return atom.dpos[1]
    def dz(atom): return atom.dpos[2]
    return eval(name)

def avg_plot(argv):
    options, args = parse_options(argv)
    configuration = process_input(args[0], options)
    output_filename = args[1]
    yfuncs = [get_atom_func(i) for i in args[2].split(",")]
    assert args[3] == "vs"
    xfunc = get_atom_func(args[4])
    xy = []
    if len(yfuncs) == 1:
        yfunc = yfuncs[0]
        for i in configuration.atoms:
            xy.append((xfunc(i), yfunc(i)))
    else:
        for i in configuration.atoms:
            xy.append((xfunc(i),) + tuple(yfunc(i) for yfunc in yfuncs))
    minx = min(i[0] for i in xy)
    maxx = max(i[0] for i in xy) + 1e-6
    print "n=%d, x E <%g,%g)" % (len(xy), minx, maxx)
    for n, yf in enumerate(yfuncs):
        miny = min(i[n+1] for i in xy)
        maxy = max(i[n+1] for i in xy)
        avgy = sum(i[n+1] for i in xy) / len(xy)
        print "%s E <%g, %g>, avg: %g" % (yf.__name__, miny, maxy, avgy)

    # make histogram
    nbins = (int(args[5]) if len(args) >= 6 else 128)
    t = nbins / (maxx - minx)
    lines = [[minx + (n + 0.5) / t] for n in range(nbins)]

    for yn in range(len(yfuncs)):
        data = [StdDev() for i in range(nbins)]
        for i in xy:
            x, y = i[0], i[1+yn]
            bin = int((x - minx) * t)
            data[bin].add_x(y)

        for n, d in enumerate(data):
            if d.n > 0:
                lines[n] += [d.mean, d.n, d.get_stddev()]

    ofile = open_any(output_filename, 'w')
    for line in lines:
        if len(line) > 1:
            print >>ofile, " ".join(str(i) for i in line)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == "--dlpoly-history-info":
        dlpoly_history_info(file(sys.argv[2]))
    elif "vs" in sys.argv:
        avg_plot(sys.argv[1:])
    else:
        convert(sys.argv[1:])


