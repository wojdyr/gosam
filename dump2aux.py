#!/usr/bin/env python
# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2
"""\
converting LAMMPS dump files AtomEye cfg format,
calculating GB energy.
Input/output files can be gzipped or bzip2'ed.
"""

usage_string = """\
Usage: 
  dump2aux.py lammps_dump output.cfg
     convert LAMMPS dump file to cfg file

  dump2aux.py lammps_dump [[energy.aux] histogram.xy]

  dump2aux.py ey energy_vs_y.histogram 
        writes GB energy vs y to gbe_vs_y.hist
"""

import sys
import bz2
import gzip

from rotmat import StdDev
from mdprim import AtomVF

#e0 = -6.1637
e0 = -6.1646668 #SiC285.tersoff
id_pos = 0
val_pos = -1
y_pos = 3
z_pos = 4
nbins = 128

atomeye_species = { 1: "12.01\nC",
                    2: "28.09\nSi",
                    3: "20.0\nB",
                    4: "20.0\nGe"
                  }

conversion_eV_A_to_J_m2 = 16.021765


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
            id_, type, x_, y_, z_ = dr.read_atom_line().split()[:5]
            pos = (float(x_) - dr.pbc_lo[0],
                   float(y_) - dr.pbc_lo[1],
                   float(z_) - dr.pbc_lo[2])
            n = int(id_)
            name = atomeye_species[int(type)].split()[-1]
            atoms[n] = AtomVF(name, n+1, pos, vel=None, force=None)
        pbc = [[self.pbc[0],0,0], [0,self.pbc[1],0], [0,0,self.pbc[2]]]
        title = "from LAMMPS dump :" + filename
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
    pos0 = 0.
    for (type, x, y, z, ex) in alist:
        if type != prev:
            cfg.write("%s\n" % atomeye_species[int(type)])
            prev = type
        xcol = (2 * (x - pos0)) % 1.
        cfg.write("%.7f %.7f %.7f %s %.3f\n" % (x, y, z, ex, xcol))


def calculate_gb_energy(dump_filename, aux_filename=None, hist_filename=None):
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
            gb_energy = delta / area * conversion_eV_A_to_J_m2
            hist_y.append((z, gb_energy, s / count))

    # the GB is assumed to be at z=0
    qb = len(hist_y) // 4
    gbe = sum(i[1] for i in hist_y[-qb:] + hist_y[:qb])
    print "GB energy: ", round(gbe, 4)

    if aux_filename:
        aux = open_any(aux_filename, "w")
        for i in values:
            print >>aux, i

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
    assert len(sys.argv) in (2,3,4), usage_string
    if sys.argv[1] == "ey":
        calc_gbe_vs_y(sys.argv[2])
    elif (len(sys.argv) == 3 and ".cfg" in sys.argv[2]):
        dump2cfg(sys.argv[1], sys.argv[2])
    else:
        calculate_gb_energy(*sys.argv[1:])


