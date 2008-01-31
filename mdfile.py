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
"""

import sys
import operator
from optparse import OptionParser
import numpy
from numpy import linalg

from mdprim import AtomVF
from model import Model
import pse


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
    return Model(atoms, pbc=[])


def import_xmol(ifile):
    number_of_atoms = int(ifile.readline())
    title = ifile.readline().strip()
    atoms = []
    for i in ifile:
        s = i.split()
        pos = (float(s[1]), float(s[2]), float(s[3]))
        atoms.append(AtomVF(s[0], len(atoms), pos, None, None))
    return Model(atoms, pbc=[], title=title)


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
    return Model(atoms, pbc, title), next_line


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


def export_for_atomeye(configuration, f, aux=None):
    "AtomEye Extended CFG format"
    if not aux:
        aux = []
    # aux = ["temperature [K]"]
    aux_fun = {
            "temperature [K]": (lambda i: i.get_temperature())
            }
    pbc = configuration.pbc 
    if pbc is None or len(pbc) == 0:
        raise ValueError("no PBC")
    if not isinstance(pbc, numpy.ndarray):
        pbc = numpy.array(pbc)
    print >>f, "Number of particles = %i" % len(configuration.atoms) 
    print >>f, "#", configuration.title
    print >>f, "A = 1.0 Angstrom (basic length-scale)"
    for i in range(3):
        for j in range(3):
            print >>f, "H0(%i,%i) = %f A" % (i+1, j+1, configuration.pbc[i][j])
    print >>f, ".NO_VELOCITY."
    entry_count = 3 + len(aux)
    entry_f = " ".join(["%f"] * entry_count)
    print >>f, "entry_count = %i" % entry_count
    for n, aux_name in enumerate(aux):
        print >>f, "auxiliary[%i] = %s" % (n, "temperature [K]")
    H_1 = linalg.inv(pbc) 
    previous_name = None
    for i in configuration.atoms:
        if previous_name != i.name:
            print >>f, pse.get_atom_mass(i.name) 
            print >>f, i.name
            previous_name = i.name
        s = numpy.dot(i.pos, H_1)
        entries = [s[0], s[1], s[2]]
        for au in aux:
            entries += aux_fun[au](i)
        print >>f, entry_f % tuple(entries)


def import_atomeye(ifile):
    # read header
    pbc = [[0,0,0], [0,0,0], [0,0,0]]
    has_velocity = True
    for line in ifile:
        if not line or line[0] == '#':
            continue
        if line.startswith("Number of particles"):
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
        if not line or line[0] == '#':
            continue
        elif line[0].isalpha():
            spec = line.strip()
        elif len(line) > 20:
            s = line.split()
            pos = (float(s[0]), float(s[1]), float(s[2]))
            pos = numpy.dot(pos, H)
            if has_velocity:
                vel = (float(s[3]), float(s[4]), float(s[5]))
                # if velocities are also reduced:
                vel = numpy.dot(vel, H) 
            atoms.append(AtomVF(spec, len(atoms), pos, vel, None))
    return Model(atoms, pbc=pbc, title="from cfg")


def import_lammps_data(ifile):
    pbc = [[0,0,0], [0,0,0], [0,0,0]]
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

    return Model(atoms, pbc=pbc, title="from lammps data")


def get_type_from_filename(name):
    lname = name.lower()
    if name.endswith(".cfg"):
        return "atomeye"
    if name.endswith(".lammps"):
        return "lammps"
    elif name.endswith(".xyz"):
        return "xmol"
    elif name.endswith(".at"):
        return "pielaszek"
    #TODO strip directory(?)
    elif "revcon" in lname or "config" in lname:
        return "dlpoly"
    elif "history" in lname:
        return "dlpoly_history"
    else:
        print "Can't deduce filetype from filename:", name
        return None


def import_autodetected(filename):
    input_type = get_type_from_filename(filename)
    infile = file(filename)
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
    else:
        assert 0


def convert():
    "converts atomistic files" 
    parser = OptionParser("usage: %prog [--pbc=list] input_file output_file")
    parser.add_option("--pbc", help="PBC, eg. '[(60,0,0),(0,60,0),(0,0,60)]'")
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("Two argumenets (input and output filenames) are required")
    input_filename = args[0]
    output_filename = args[1]
    configuration = import_autodetected(input_filename)
    if options.pbc:
        configuration.pbc = eval(options.pbc)
    output_type = get_type_from_filename(output_filename)
    #TODO: check if file already exists
    ofile = file(output_filename, 'w')
    configuration.export_atoms(output_filename, output_type)


if __name__ == '__main__':
    if sys.argv[1] == "--dlpoly-history-info":
        dlpoly_history_info(file(sys.argv[2]))
    else:
        convert()



