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
           * miller indices prefixed with letter m (e.g. m011) meaning
             median plane; the boundary will be calculated as the median plane
             rotated by theta/2 around the axis.

 - instead of sigma (one number) you can give:
   * m,n (e.g. 23,4)
   * theta=value, i.e. value of angle in degrees (e.g. theta=90)
 - dim_x, dim_y and dim_z are in nm
 - options:
   * nofit - if this option is _not_ specified, PBC dimensions will be tuned
             to make the system periodic
   * mono1 - generate only upper half of the bicrystal, i.e. monocrystal
   * mono2 - generate only bottom half of the bicrystal, i.e. monocrystal
   * remove:dist - If there are two atoms in distance < dist [Angstroms],
             one of the atoms is removed.
   * remove2:dist - for binary systems only; like the option above, but only
             pairs of atoms of the same species are checked.
   * vacuum:length - vacuum in z direction. Makes 2D slab with z dimension
             increased by the length.
   * shift:dx,dy,dz - shift nodes in unit cell.
   * lattice:name - e.g. SiC. For heterophase boundaries, two names should be
                    separated by comma without spaces (e.g. lattice:Cu,Fe).
                    In this case one of the grains may not fit well into PBC.
   * edge:z1,z2 - Removes atoms that have y in lower half of the box and 
             z1 < z < z2. There is a chance that this will become an edge 
             dislocation after squeezing or running high temperature MD.

Examples:
    bicrystal.py 001 twist 5 20 20 80 twist_s5.cfg
    bicrystal.py 100 m011 5 20 20 80 tilt_s5.cfg
    bicrystal.py 100 0,1,0 theta=90 1 1 1  tmp.cfg

See the wiki for more examples.

caution: the program was tested only for a few cases (may not work in others)
"""

import math
from math import sqrt, degrees, radians, acos, ceil
import sys
from copy import deepcopy
import numpy
from numpy import dot, identity, inner, zeros
from numpy import linalg
from monocryst import RotatedMonocrystal, OrthorhombicPbcModel, \
                     get_command_line, get_named_lattice
import csl
from rotmat import rodrigues, print_matrix, round_to_multiplicity


class Bicrystal(OrthorhombicPbcModel):
    def __init__(self, lattice1, lattice2, dim, rot_u, rot_b, title):
        OrthorhombicPbcModel.__init__(self, lattice1, dim, title=title)
        self.mono_u = RotatedMonocrystal(lattice1, dim, rot_u)
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
        self.vacuum = None # margin for dim z
        self.dim = None
        self.fit = None
        self.zfit = None
        self.mono1 = None
        self.mono2 = None
        self.remove_dist = None
        self.remove_dist2 = None
        self.all = None
        self.allall = None
        self.lattice_name = "sic"
        self.lattice_shift = None
        self.edge = None
        self.antiphase = False # unused


    def parse_sigma_and_find_theta(self, sigma_arg):
        if sigma_arg.startswith("theta="):
            sigma = None
            m, n = None, None
            theta = radians(float(sigma_arg[6:]))
        elif "," not in sigma_arg:
            if sigma_arg.startswith("u"):
                sigma = int(sigma_arg[1:])
                min_angle = radians(45.)
            else:
                sigma = int(sigma_arg)
                min_angle = None
            r = csl.find_theta(self.axis, sigma, min_angle=min_angle)
            if r is None:
                print "CSL not found! Wrong sigma or axis?"
                sys.exit()
            theta, m, n = r
        else:
            m_, n_ = sigma_arg.split(",")
            m, n = int(m_), int(n_)
            sigma = csl.get_cubic_sigma(self.axis, m, n)
            theta = csl.get_cubic_theta(self.axis, m, n)
        if sigma is not None:
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
        fit_dim = []
        if self.fit:
            fit_dim += [0, 1]
            if self.zfit:
                fit_dim += [2]
        for i in fit_dim:
            mult = ceil(float(dim[i]) / min_dim[i]) or 1
            dim[i] = mult * min_dim[i]
            # dim[i] = round_to_multiplicity(min_dim[i], dim[i])
        if self.vacuum:
            dim[2] += self.vacuum # margin in dim z
        print "-------> dimensions [A]: ", dim[0], dim[1], dim[2]
        self.dim = dim


def print_boundary_type(axis, plane, theta):
    if (plane == axis * sum(plane) / sum(axis)).all():
        bt = "twist"
    elif inner(plane, axis) == 0:
        R = rodrigues(axis, theta, verbose=False)
        p2 = dot(linalg.inv(R), plane)
        plane2 = csl.scale_to_integers(p2)
        bt = "tilt (%i %i %i) (%i %i %i)" % (tuple(plane) + tuple(plane2))
    else:
        bt = "mixed"
    return bt


def parse_args():
    if len(sys.argv) < 7:
        print usage
        sys.exit()

    opts = BicrystalOptions()
    opts.axis = csl.parse_miller(sys.argv[1])
    print "-------> rotation axis: [%i %i %i]" % tuple(opts.axis)

    opts.parse_sigma_and_find_theta(sys.argv[3])

    plane = sys.argv[2]
    if plane == "twist":
        opts.plane = opts.axis.copy()
    elif plane.startswith("m"):
        m_plane = csl.parse_miller(plane[1:])
        if inner(m_plane, opts.axis) != 0:
            raise ValueError("Axis must be contained in median plane.")
        R = rodrigues(opts.axis, opts.theta / 2., verbose=False)
        plane_ = dot(R, m_plane)
        opts.plane = csl.scale_to_integers(plane_)
    else:
        opts.plane = csl.parse_miller(plane)
    print "-------> boundary plane: (%i %i %i)" % tuple(opts.plane)
    bt = print_boundary_type(opts.axis, opts.plane, opts.theta)
    print "         boundary type:", bt

    opts.req_dim = [float(eval(i, math.__dict__)) for i in sys.argv[4:7]]

    options = sys.argv[7:-1]
    for i in options:
        if i == "nofit":
            assert opts.fit is None
            opts.fit = False
        if i == "nozfit":
            assert opts.zfit is None
            opts.zfit = False
        elif i == "mono1":
            assert opts.mono1 is None
            opts.mono1 = True
        elif i == "mono2":
            assert opts.mono2 is None
            opts.mono2 = True
        elif i == "all":
            assert opts.all is None and opts.allall is None
            opts.all = True
        elif i == "allall":
            assert opts.allall is None and opts.all is None
            opts.allall = True
        elif i.startswith("remove:"):
            assert opts.remove_dist is None
            opts.remove_dist = float(i[7:])
        elif i.startswith("remove2:"):
            assert opts.remove_dist2 is None
            opts.remove_dist2 = float(i[8:])
        elif i.startswith("vacuum:"):
            assert opts.vacuum is None
            opts.vacuum = float(i[7:]) * 10. #nm -> A
        elif i.startswith("lattice:"):
            opts.lattice_name = i[8:]
        elif i.startswith("shift:"):
            s = i[6:].split(",")
            if len(s) != 3:
                raise ValueError("Wrong format of shift parameter")
            opts.lattice_shift = [float(i) for i in s]
        elif i.startswith("edge:"):
            try:
                z1, z2 = i[5:].split(",")
                opts.edge = (float(z1), float(z2))
            except (TypeError, ValueError):
                raise ValueError("Wrong format of edge parameter")
        else:
            raise ValueError("Unknown option: %s" % i)
    # default values
    if opts.fit is None:
        opts.fit = True
    if opts.zfit is None:
        opts.zfit = True
    if opts.mono1 is None:
        opts.mono1 = False
    if opts.mono2 is None:
        opts.mono2 = False
    #if opts.remove_dist is None:
    #    opts.remove_dist = 0.8 * opts.atom_min_dist

    opts.output_filename = sys.argv[-1]
    return opts


def main():
    opts = parse_args()

    # R is a matrix that transforms lattice in the bottom monocrystal
    # to lattice in the upper monocrystal
    R = rodrigues(opts.axis, opts.theta)

    if opts.sigma:
        # C is CSL primitive cell
        C = csl.find_csl_matrix(opts.sigma, R)
        print_matrix("CSL primitive cell", C)

        ## and now we determine CSL for fcc lattice
        #C = csl.pc2fcc(C)
        #C = csl.beautify_matrix(C)
        #print_matrix("CSL cell for fcc:", C)
    else:
        C = identity(3)

    # CSL-lattice must be periodic is our system.
    # * PBC box must be orthonormal
    # * boundaries must be perpendicular to z axis of PBC box

    print "CSL cell with z || [%s %s %s]" % tuple(opts.plane),
    Cp = csl.make_parallel_to_axis(C, col=2, axis=opts.plane)
    print_matrix("", Cp)

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
    assert (numpy.abs(invrot - linalg.inv(rot)) < 1e-9).all(), "%s != %s" % (
                                                     invrot, linalg.inv(rot))
    #print "hack warning: min_dim[1] /= 2."
    #min_dim[1] /= 2.

    if "," in opts.lattice_name:
        name1, name2 = opts.lattice_name.split(",")
        lattice1 = get_named_lattice(name1)
        lattice2 = get_named_lattice(name2)
    else:
        lattice1 = get_named_lattice(opts.lattice_name)
        lattice2 = deepcopy(lattice1)

    if opts.lattice_shift:
        lattice1.shift_nodes(opts.lattice_shift)
        lattice2.shift_nodes(opts.lattice_shift)

    # the anti-phase GB can be made by swapping two species in one grain
    if opts.antiphase:
        assert lattice2.count_species() == 2
        lattice2.swap_node_atoms_names()

    a = lattice1.unit_cell.a
    opts.find_dim([i * a for i in min_dim])

    #rot_mat1 = rodrigues(opts.axis, rot1)
    #rot_mat2 = rodrigues(opts.axis, rot2)
    rot_mat1 = dot(linalg.inv(R), invrot)
    rot_mat2 = invrot
    #print "rot1", "det=%g" % linalg.det(rot_mat1), rot_mat1
    #print "rot2", "det=%g" % linalg.det(rot_mat2), rot_mat2

    title = get_command_line()
    if opts.mono1:
        config = RotatedMonocrystal(lattice1, opts.dim, rot_mat1,
                                    title=title)
    elif opts.mono2:
        config = RotatedMonocrystal(lattice2, opts.dim, rot_mat2,
                                    title=title)
    else:
        config = Bicrystal(lattice1, lattice2, opts.dim, rot_mat1, rot_mat2,
                           title=title)
    config.generate_atoms(z_margin=opts.vacuum)

    if not opts.mono1 and not opts.mono2 and opts.remove_dist > 0:
        print "Removing atoms in distance < %s ..." % opts.remove_dist
        config.remove_close_neighbours(opts.remove_dist)

    if opts.remove_dist2:
        a_atoms = []
        b_atoms = []
        a_name = config.atoms[0].name
        for i in config.atoms:
            if i.name == a_name:
                a_atoms.append(i)
            else:
                b_atoms.append(i)
        config.atoms = []
        for aa in a_atoms, b_atoms:
            print "Removing atoms where %s-%s distance is < %s ..." % (
                                     aa[0].name, aa[0].name, opts.remove_dist2)
            config.remove_close_neighbours(distance=opts.remove_dist2, atoms=aa)
            config.atoms += aa

    if opts.edge:
        z1, z2 = opts.edge
        sel = [n for n, a in enumerate(config.atoms)
               if 0 < a.pos[1] <= config.pbc[1][1] / 2. and z1 < a.pos[2] < z2]
        print "edge: %d atoms is removed" % len(sel)
        for i in reversed(sel):
            del config.atoms[i]

    if opts.all:
        #config.output_all_removal_possibilities(opts.output_filename)
        config.apply_all_possible_cutoffs_to_stgb(opts.output_filename,
                                                  single_cutoff=True)
        return

    if opts.allall:
        config.apply_all_possible_cutoffs_to_stgb(opts.output_filename,
                                                  single_cutoff=False)
        return

    config.export_atoms(opts.output_filename)


if __name__ == '__main__':
    main()


