#!/usr/bin/env python
# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2
"""\
Converts LAMMPS dump files to AtomEye cfg format,
calculates GB energy in bicrystal geometry.
Input/output files can be gzipped or bzip2'ed.
"""

usage_string = """\
Usage:
  dump2aux.py lammps_dump output.cfg
     convert LAMMPS dump file to cfg file

  dump2aux.py hist lammps_dump histogram.xy
        this may not work now

  dump2aux.py ey energy_vs_y.histogram
        write GB energy vs y to gbe_vs_y.hist

  dump2aux.py lammps_dump1 [lammps_dump2 ...]
        calculate GB energies
"""

import sys
import os.path
import bz2
import gzip
from math import sqrt

from rotmat import StdDev
from mdprim import AtomVF
import model

#e0 = -6.1637
e0 = -6.1646668 #SiC285.tersoff
id_pos = 0
val_pos = -1
y_pos = 3
z_pos = 4
nbins = 128
#gb_relative_width = 0.7
gb_relative_width = None

atomeye_species = { 1: "12.01\nC",
                    2: "28.09\nSi",
                    3: "20.0\nB",
                    4: "20.0\nCl"
                  }

conversion_eV_A2_to_J_m2 = 16.021765


def open_any(name, mode='r'):
    if name.endswith(".bz2"):
        return bz2.BZ2File(name, mode)
    elif name.endswith(".gz"):
        return gzip.GzipFile(name, mode)
    else:
        return file(name, mode)


class DumpReader:
    def __init__(self, dump_filename):
        self.dump = open_any(dump_filename)
        self.filename = dump_filename
        self._read_header()

    def _read_header(self):
        assert self.dump.readline().strip() == "ITEM: TIMESTEP"
        self.dump.readline() # skip timestep
        assert self.dump.readline().strip() == "ITEM: NUMBER OF ATOMS"
        self.natoms = int(self.dump.readline())
        assert self.dump.readline().strip() == "ITEM: BOX BOUNDS"
        self.pbc = []
        self.pbc_lo = []
        for i in range(3):
            pmin_, pmax_ = self.dump.readline().split()
            pmin, pmax = float(pmin_), float(pmax_)
            self.pbc.append(pmax - pmin)
            self.pbc_lo.append(pmin)
        atoms_line = self.dump.readline()
        assert atoms_line.startswith("ITEM: ATOMS id type x y z")
        self.extra_data = atoms_line.split()[7:]

    def read_atom_line(self):
        return self.dump.readline()

    def get_configuration(self):
        atoms = [None for i in range(self.natoms)]
        for i in range(self.natoms):
            id_, type, x_, y_, z_ = self.read_atom_line().split()[:5]
            pos = (float(x_) - self.pbc_lo[0],
                   float(y_) - self.pbc_lo[1],
                   float(z_) - self.pbc_lo[2])
            n = int(id_)
            name = atomeye_species[int(type)].split()[-1]
            atoms[n-1] = AtomVF(name, n, pos, vel=None, force=None)
        pbc = [[self.pbc[0],0,0], [0,self.pbc[1],0], [0,0,self.pbc[2]]]
        title = "from LAMMPS dump :" + self.filename
        return model.Model(atoms, pbc=pbc, title=title)


def dump2cfg(dump_filename, cfg_filename):
    "converts LAMMPS dump_filename to AtomEye extended CFG file"
    # read header
    dr = DumpReader(dump_filename)
    #natoms, pbc, pbc_lo = read_dump_header(dump)
    # write header
    cfg = open_any(cfg_filename, "w")
    cfg.write("Number of particles = %d\n" % dr.natoms);
    cfg.write("# converted by dump2aux.dump2cfg from %s\n" % dump_filename);
    cfg.write("# full original path: %s\n" % os.path.abspath(dump_filename));
    cfg.write("A = 1.0 Angstrom (basic length-scale)\n")
    for i in range(3):
        for j in range(3):
            cfg.write("H0(%i,%i) = %f A\n" % (i+1, j+1,
                                              (dr.pbc[i] if i == j else 0.)))
    cfg.write(".NO_VELOCITY.\n")
    cfg.write("entry_count = %d\n" % (4+len(dr.extra_data)))
    for n, name in enumerate(dr.extra_data):
        cfg.write("auxiliary[%d] = %s\n" % (n, name))
    cfg.write("auxiliary[%d] = xcolor [0-1]\n" % len(dr.extra_data))
    # read atoms
    alist = [None for i in range(dr.natoms)]
    for i in range(dr.natoms):
        id_, type, x_, y_, z_, ex = dr.read_atom_line().split(None, 5)
        x = (float(x_) - dr.pbc_lo[0]) / dr.pbc[0]
        y = (float(y_) - dr.pbc_lo[1]) / dr.pbc[1]
        z = (float(z_) - dr.pbc_lo[2]) / dr.pbc[2]
        alist[int(id_)-1] = (type, x, y, z, ex.strip())
    # write atoms
    prev = None
    pos0 = _find_pos0(alist)
    for (type, x, y, z, ex) in alist:
        if type != prev:
            cfg.write("%s\n" % atomeye_species[int(type)])
            prev = type
        xcol = (2 * (x - pos0)) % 1.
        cfg.write("%.7f %.7f %.7f %s %.3f\n" % (x, y, z, ex, xcol))

# pos0 is used for xcolor, i.e. for coloring based on the x coordinate
# selection of pos0 matters when we compare different frames
# If there is a rigid surface in the system, we use it as a reference for
# other coordinates.
def _find_pos0(alist):
    selected = [a for a in alist if a[0] == "3"]
    if selected:
        elem = min(selected, key = lambda x: x[3])
        return elem[1]
    else: # no rigid surface
        return sum(a[1] for a in alist) / len(alist)

def _print_gb_energy(energies, dr, verbose=False):
    count = len(energies)
    energy = sum(energies)
    excess = energy - count * e0
    if verbose:
        print "PBC: %g x %g x %g" % tuple(dr.pbc)
        print "total energy of %d atoms: %g eV (%g/at.), excess: %g" % (
                count, energy, energy/count, excess)
    area = dr.pbc[0] * dr.pbc[1]
    gb_energy = excess / area * conversion_eV_A2_to_J_m2
    if verbose:
        print "GB energy:", gb_energy
        print "n*Edisl:", excess * conversion_eV_A2_to_J_m2 * 1e-10 / dr.pbc[0]
    return gb_energy

def calculate_dislocation_energy(dump_filename, y0, z0, r):
    dr = DumpReader(dump_filename)
    energies = []
    for i in range(dr.natoms):
        s = dr.read_atom_line().split()
        dy = abs(float(s[3]) - y0)
        if dy > dr.pbc[1] / 2:
            dy = dr.pbc[1] - dy
        dz = abs(float(s[4]) - z0)
        if dz > dr.pbc[2] / 2:
            dz = dr.pbc[2] - dz
        if dy*dy + dz*dz < r*r:
            energies.append(float(s[val_pos]))
    count = len(energies)
    energy = sum(energies)
    excess = energy - count * e0
    disl_length = dr.pbc[0]
    eV_A_to_J_m = conversion_eV_A2_to_J_m2 * 1e-10
    e_disl = excess / disl_length * eV_A_to_J_m
    print "%d atoms in cyllinder (y-%g)^2 + (z-%g)^2 < %g^2" % (count,y0,z0,r)
    print "Edisl [J/m]: %g" % (e_disl)
    return e_disl

def calculate_gbe_of_types12(dump_filename):
    # read and parse first snapshot
    dr = DumpReader(dump_filename)

    energies = []
    for i in range(dr.natoms):
        tokens = dr.read_atom_line().split()
        if tokens[1] in ("1", "2"):
            val = float(tokens[val_pos])
            energies.append(val)

    return _print_gb_energy(energies, dr)


def calculate_total_energy(dump_filename):
    dr = DumpReader(dump_filename)

    energies = []
    for i in dr.dump:
        s = i.split()
        y = float(s[3]) / dr.pbc[1] % 1.
        z = float(s[4]) / dr.pbc[2] % 1.
        #if not 0.25 <= y < 0.75:
        #if not 0.30 <= z < 0.70:
        if 1:
            energies.append(float(s[val_pos]))
    #energies = [float(i.split()[val_pos]) for i in dr.dump]

    _print_gb_energy(energies, dr, verbose=True)


def calculate_gb_energy(dump_filename, hist_filename=None):
    if gb_relative_width is None:
        return calculate_gbe_of_types12(dump_filename)
        ########## that's the end ####################

    # read and parse first snapshot
    dr = DumpReader(dump_filename)
    area = dr.pbc[0] * dr.pbc[1]

    values = [None for i in range(dr.natoms)]
    hist = [[] for i in range(nbins)]

    for i in range(dr.natoms):
        tokens = dr.read_atom_line().split()
        id = int(tokens[id_pos])
        val = float(tokens[val_pos])
        values[id-1] = val

        z = float(tokens[z_pos]) % dr.pbc[2]
        bin = int(nbins * z / dr.pbc[2])
        assert bin < nbins
        hist[bin].append(values[id-1])

    #print "total energy:", sum(values)

    hist_y = []
    for n, vv in enumerate(hist):
        if vv:
            z = (n + 0.5) / nbins * dr.pbc[2]
            s = sum(vv)
            count = len(vv)
            delta = s - e0 * count
            gb_energy = delta / area * conversion_eV_A2_to_J_m2
            hist_y.append((z, gb_energy, s / count))

    # the GB is assumed to be at z=0
    qb = int(gb_relative_width * len(hist_y) / 2)
    gbe = sum(i[1] for i in hist_y[-qb:] + hist_y[:qb])

    #if aux_filename:
    #    aux = open_any(aux_filename, "w")
    #    for i in values:
    #        print >>aux, i

    if hist_filename:
        hist_file = open_any(hist_filename, "w")
        for (z, gb_energy, avg) in hist_y:
            print >>hist_file, z, avg, gb_energy

    return gbe


def calc_gbe_vs_y(dump_filename):
    "writes GB energy vs y histogram to file gbe_vs_y.hist"
    # read and parse first snapshot
    nbins = 50
    dr = DumpReader(dump_filename)
    width = 10.0
    Y = dr.pbc[1]
    area = dr.pbc[0] * dr.pbc[1]

    hist = [0.0 for i in range(nbins)]
    for i in range(dr.natoms):
        tokens = dr.read_atom_line().split()
        z = float(tokens[z_pos])
        if abs(z) > width and abs(z) < dr.pbc[2] - width:
            continue
        y = float(tokens[y_pos]) % Y
        bin = int(nbins * y / Y)
        assert bin < nbins
        delta = float(tokens[val_pos]) - e0
        hist[bin] += delta
    print "GBE:", sum(hist) / area * conversion_eV_A_to_J_m2

    hist_file = open_any("gbe_vs_y.hist", "w")
    for n, d in enumerate(hist):
        print >>hist_file, (n+0.5) / nbins, d / area * conversion_eV_A_to_J_m2


if __name__ == "__main__":
    assert len(sys.argv) >1, usage_string
    if len(sys.argv) == 3 and sys.argv[1] == "ey":
        calc_gbe_vs_y(sys.argv[2])
    if len(sys.argv) >= 3 and sys.argv[1] == "total":
        for f in sys.argv[2:]:
            if len(sys.argv) > 3:
                print f
            calculate_total_energy(f)
    elif len(sys.argv) == 4 and sys.argv[1] == "hist":
        gbe = calculate_gb_energy(sys.argv[2], sys.argv[3])
        print "GB energy: ", round(gbe, 4)
    elif len(sys.argv) == 3 and ".cfg" in sys.argv[2]:
        dump2cfg(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 6 and sys.argv[1] == "disl":
        calculate_dislocation_energy(dump_filename=sys.argv[2],
                                     y0=float(sys.argv[3]),
                                     z0=float(sys.argv[4]),
                                     r=float(sys.argv[5]))
    else:
        if gb_relative_width:
            print "GB energy [J/m2]"
            print "GB width assumed as %g%% of slab" % (gb_relative_width * 100)
        for i in sys.argv[1:]:
            gbe = calculate_gb_energy(i)
            name = (i[:-8] if i.endswith(".dump.gz") else i)
            print "%s\t%f" % (name, gbe)


