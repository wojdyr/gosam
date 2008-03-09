#!/usr/bin/env python
# this file is part of gosam (generator of simple atomistic models) 
# Licence: GNU General Public License version 2
"""\
tool for generating bicrystals
"""

usage = """\
Usage: 
   bicrystal.py axis plane sigma dim_x dim_y dim_z [options] output_file

 - axis of rotation should be given as three numbers, e.g.: "001", "111"
 - boundary is always at plane z = (dim_z / 2)
 - plane - it can be given as:
           * miller indices of the boundary plane in bottom monocrystal lattice
           * "twist" - keyword that means that plane is perpendicular to axis
           * TODO: miller indices prefixed with letter s (e.g. s001) meaning
             symmetry place (the boundary will be calculated)
             
 - instead of sigma (one number) you can give m,n (e.g. 23,4)
 - dim_x, dim_y and dim_z are in nm
 - options:
   * nofit - if this option is _not_ specified, PBC dimensions will be tuned 
             to make the system periodic 
   * mono1 - generate only upper half of the bicrystal, i.e. monocrystal
   * mono2 - generate only bottom half of the bicrystal, i.e. monocrystal
   * remove:dist - by default, if there are two atoms in distance < 80% of
             minimal interatomic distance in perfect crystal, one of the atoms
             is removed. You can specify "remove:0" to disable atom removing,
             or other value [in Angstroms] to change the default distance.

Examples:
    bicrystal.py 001 twist 5 20 20 80 twist_s5.cfg 
    bicrystal.py 100 013 5 20 20 80 tilt_s5.cfg 
    TODO: bicrystal.py 100 s001 5 20 20 80 tilt_s5.cfg 

caution: the program was tested only for a few cases (may not work in others)

"""


import math
from math import sin, cos, pi, atan, sqrt, degrees, radians, asin, acos
import sys
from copy import deepcopy
import random
import numpy
from numpy import dot, array, identity, inner, zeros
from numpy import linalg
from monocryst import RotatedMonocrystal, OrthorhombicPbcModel, \
                     make_sic_lattice, get_command_line
import csl
from rotmat import rodrigues, print_matrix, round_to_multiplicity


random.seed(12345)


class Bicrystal(OrthorhombicPbcModel):
    def __init__(self, lattice, dim, rot_u, rot_b, title):
        OrthorhombicPbcModel.__init__(self, lattice, dim, title=title)
        self.mono_u = RotatedMonocrystal(deepcopy(lattice), dim, rot_u)
        lattice2 = deepcopy(lattice)

        # we want stechiometry in binary systems, and this helps:
        # update: I'm not sure if the bicrystal will be still "CSL" after
        # this trick, and if each boundary (not only total) will be 
        # stechiometric
        #if lattice2.count_species() == 2:
        #    lattice2.swap_node_atoms_names()

        self.mono_b = RotatedMonocrystal(lattice2, dim, rot_b)


    def generate_atoms(self, z_margin=0.):
        #print "Bicrystal.generate_atoms"
        self.atoms = (self.mono_u.generate_atoms(upper=True, z_margin=z_margin)
                  + self.mono_b.generate_atoms(upper=False, z_margin=z_margin))
        print "Number of atoms in bicrystal: %i" % len(self.atoms) 
        self.print_boundary_angle()

    def print_boundary_angle(self):
        def calc_angle(v1, v2):
            return acos( dot(v1, v2) / sqrt(inner(v1,v1) * inner(v2,v2)) )
        u0 = self.mono_u.unit_cell.get_unit_shift(0)
        b = [self.mono_b.unit_cell.get_unit_shift(i) for i in range(3)]
        b += [-i for i in b]
        angles = [degrees(calc_angle(u0, i)) for i in b]
        print "angles between upper and bottom:", \
                ", ".join("%.1f" % i for i in angles)


class BicrystalOptions:
    def __init__(self):
        self.axis = None
        self.plane = None
        self.sigma = None
        self.theta = None
        self.m = None
        self.n = None
        self.req_dim = None
        self.margin = 0.0 # margin for dim z
        self.dim = None
        self.fit = None
        self.mono1 = None
        self.mono2 = None
        self.remove_dist = None

        self.lattice = make_sic_lattice()
        a = self.lattice.unit_cell.a
        # for zinc blende structure:
        self.atom_min_dist = a * sqrt(3) / 4 


    def parse_sigma_and_find_theta(self, sigma_arg):
        if "," not in sigma_arg:
            sigma = int(sigma_arg)
            r = csl.find_theta(self.axis, sigma)
            if r is None:
                print "CSL not found! Wrong sigma or axis?"
                sys.exit()
            theta, m, n = r
        else:
            m_, n_ = sigma_arg.split(",")
            m, n = int(m_), int(n_)
            sigma, theta = csl.cubic_csl(self.axis, m, n)
        print "-------> sigma = %i" % sigma
        print "-------> theta = %.3f deg" % degrees(theta)
        self.sigma = sigma
        self.theta = theta
        self.m = m
        self.n = n


    def find_dim(self, min_dim):
        print "-------> min. dim.  [A]: ", min_dim[0], min_dim[1], min_dim[2]
        dim = [i * 10 for i in self.req_dim] # nm -> A
        if self.mono1 or self.mono2:
            dim[2] *= 2
        dim[2] += 2 * self.margin # margin in dim z
        if self.fit:
            #min_dim = self.find_min_dim()
            dim = [round_to_multiplicity(min_dim[n], i) 
                                                   for n, i in enumerate(dim)]
        print "-------> dimensions [A]: ", dim[0], dim[1], dim[2]
        self.dim = dim


def parse_miller(s):
    if len(s) == 3 and s.isdigit():
        return array([int(s[i]) for i in range(3)])
    elif ',' in s:
        sp = s.split(",")
        assert len(sp) == 3
        return array([int(i) for i in sp])
    else:
        raise ValueError("Can't parse miller indices: %s" % s)


def parse_args():
    if len(sys.argv) < 7:
        print usage
        sys.exit()

    opts = BicrystalOptions()
    opts.axis = parse_miller(sys.argv[1])
    print "-------> rotation axis: [%i %i %i]" % tuple(opts.axis)

    plane = sys.argv[2]
    if plane == "twist":
        opts.plane = opts.axis.copy()
    elif plane.startswith("s"):
        s_plane = parse_miller(plane[1:])
        assert 0, "TODO"
        #opts.plane = 
    else:
        opts.plane = parse_miller(plane)
    print "-------> boundary plane: (%i %i %i)" % tuple(opts.plane)

    opts.parse_sigma_and_find_theta(sys.argv[3])

    opts.req_dim = [float(eval(i, math.__dict__)) for i in sys.argv[4:7]]
    options = sys.argv[7:-1]
    opts.fit = True
    opts.mono1 = opts.mono2 = False
    opts.remove_dist = 0.8 * opts.atom_min_dist
    for i in options:
        if i == "nofit":
            opts.fit = False
        elif i == "mono1":
            opts.mono1 = True
        elif i == "mono2":
            opts.mono2 = True
        elif i.startswith("remove:"):
            opts.remove_dist = float(i[7:])
        else:
            raise ValueError("Unknown option: %s" % i)
    opts.output_filename = sys.argv[-1]
    return opts


def main():
    opts = parse_args()

    # R is a matrix that transforms lattice in the bottom monocrystal
    # to lattice in the upper monocrystal 
    R = rodrigues(opts.axis, opts.theta)

    # C is CSL primitive cell
    C = csl.find_csl_matrix(opts.sigma, R)
    print_matrix("CSL primitive cell", C)
    
    # and now we determine CSL for fcc lattice
    C = csl.pc2fcc(C)
    C = csl.beautify_matrix(C)
    print_matrix("CSL cell for fcc:", C)

    # CSL-lattice must be periodic is our system.
    # * PBC box must be orthonormal 
    # * boundaries must be perpendicular to z axis of PBC box

    Cp = csl.make_parallel_to_axis(C, col=2, axis=opts.plane)
    print_matrix("CSL cell with z || [%s %s %s]" % tuple(opts.plane), Cp)
    Cp = C

    min_pbc = csl.find_orthorhombic_pbc(Cp)
    print_matrix("Minimal(?) orthorhombic PBC", min_pbc)

    min_dim = []
    pbct = min_pbc.transpose().astype(float)
    rot = zeros((3, 3))
    for i in range(3):
        length = sqrt(inner(pbct[i], pbct[i]))
        rot[i] = pbct[i] / length
        min_dim.append(length)
    invrot = rot.transpose()
    assert (numpy.abs(invrot - linalg.inv(rot)) < 1e-9).all()
    a = opts.lattice.unit_cell.a
    opts.find_dim([i * a for i in min_dim])

    #rot_mat1 = rodrigues(opts.axis, rot1) 
    #rot_mat2 = rodrigues(opts.axis, rot2) 
    rot_mat1 = dot(linalg.inv(R), invrot)
    rot_mat2 = invrot
    #print "rot1", rot_mat1 
    #print "rot2", rot_mat2 

    title = get_command_line()
    if opts.mono1:
        config = RotatedMonocrystal(opts.lattice, opts.dim, rot_mat1,
                                    title=title)
    elif opts.mono2:
        config = RotatedMonocrystal(opts.lattice, opts.dim, rot_mat2,
                                    title=title)
    else:
        config = Bicrystal(opts.lattice, opts.dim, rot_mat1, rot_mat2,
                           title=title)
    config.generate_atoms(z_margin=opts.margin)
    if not opts.mono1 and not opts.mono2 and opts.remove_dist > 0:
        config.remove_close_neighbours(opts.remove_dist)
    config.export_atoms(opts.output_filename)


#def get_random_rotation():
#    v = [1, 1, 1]
#    while inner(v, v) > 1:
#        for i in range(3):
#            v[i] = random.random()
#    theta = random.uniform(0, pi)
#    rot_mat = rodrigues(v, theta) 
#    return rot_mat
#
#
## outdated
#def random_bicrystal(): 
#    lattice = make_default_lattice()
#    dim = [100, 100, 200]
#    rot1 = get_random_rotation()
#    rot2 = get_random_rotation()
#    config = Bicrystal(lattice, dim, rot1, rot2)
#    config.generate_atoms()
#    cryst_min_dist = lattice.unit_cell.a * sqrt(3) / 4 # in zinc blende
#    config.remove_close_neighbours(0.8 * cryst_min_dist)
#    #config.export_atoms("bicr.xyz", format="xmol")
#    config.export_atoms("bicr.cfg", format="atomeye")



if __name__ == '__main__':
    main()


