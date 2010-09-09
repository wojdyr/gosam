#!/usr/bin/env python
# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2
# $Id$
"""
Module for geometrical modelling of crystal grain.

One of possible uses is as follows:

1  Define unit cell (eg. for triclinic - a, b, c, alpha, beta, gamma)

2  Compute parameters of unit cell volume, reciprocal lattice parameters,
    transformation matrix M (triclinic -> orthonormal axes), ...

3  Define atom positions in unit cell (nodes in cell * atoms in node)

4  Define grain surfaces and pressures - set of planes,
    each defined by (hkl) and r - distance from origin,
    and, optionally, distance D_shell and function emulating pressure
     (changing position) if distance from atom to plane is smaller then D_shell

5 Compute parameters of normal plane equations of surfaces.

6  Get min/max value along a,b and c axes
     this step includes computing halfspace intersection about a point
     (using qhull program) - it should give all vertices.

7  (optional) Show shape of grain using Geomview.

8  Iterate over all possible cells (using min/max values) and atoms.
   Checking if _node_ is inside of grain; if atoms are in grain shell
   (distance from surface < D_shell), moving atoms in shell according
   to given function (in direction normal to surface).

9 (optional) Make vacancies

9  (optional) Modify atom position (eg. to simulate thermal vibation)

10  Write output
"""

from math import cos, sin, acos, asin, sqrt, pi, floor, ceil, radians, degrees
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
import commands
import sys
import os.path
from numpy import array, zeros, dot, sometrue, inner, cross
# file formats
from mdfile import *
from latt import *
import mdprim
from model import Model
from rotmat import rodrigues


fcc_nodes = [
    (0.0, 0.0, 0.0),
    (0.5, 0.5, 0.0),
    (0.0, 0.5, 0.5),
    (0.5, 0.0, 0.5),
]

bcc_nodes = [
    (0.0, 0.0, 0.0),
    (0.5, 0.5, 0.5)
]


class UnexpectedArgsError(Exception):
    def __init__(self, value):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return repr(self.value)


class NotInitializedError(Exception):
    def __init__(self, value):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return repr(self.value)



class Plane:
    "2D plane in 3D space. Defined using 4 paramaters."
    def __init__(self, ABCD=None):
        self.initialized = False
        if ABCD:
            self.set_ABCD(ABCD)

    def __str__(self):
        return self.describe_ABCD() + "\nor, in other words,\n" \
                + self.describe_angles()


    def set_ABCD(self, ABCD):
        "set as (A,B,C,D) parameters - A*x+B*y+C*z+D=0"
        if len(ABCD) != 4:
            raise UnexpectedArgsError("4 (A B C D) arguments for plane "
                                              "required, got %s" % len(ABCD))
        self.A, self.B, self.C, self.D = ABCD
        self.initialized = True
        self._compute_angles_from_ABCD()


    def set_angles(self, angles):
        """set as (alpha,beta,gamma,p) parameters of normal plane equation:
                 x*cos(alpha)+y*cos(beta)+z*cos(gamma)-p=0"""
        if len(angles) != 4:
            raise UnexpectedArgsError("4 (alpha,beta,gamma,p) arguments"
                                    "for plane required, got %s" % len(angles))
        self.alpha, self.beta, self.gamma, self.p = angles
        self.initialized = True
        self._compute_ABCD_from_angles()


    def _compute_ABCD_from_angles(self):
        "this method is always called when parameters are changed"
        self.A,self.B,self.C = cos(self.alpha), cos(self.beta), cos(self.gamma)
        self.D = -self.p
        #some precomputations -- for optimization
        self.cosines = array((self.A, self.B, self.C))


    def _compute_angles_from_ABCD(self):
        A,B,C,D = self.A, self.B, self.C, self.D
        mi = 1./sqrt(A**2 + B**2 + C**2)
        if D>0:
            mi = -mi
        self.alpha = acos(A*mi)
        self.beta = acos(B*mi)
        self.gamma = acos(C*mi)
        self.p = -D*mi
        self._compute_ABCD_from_angles() # to ensure A,B,C,D is normalized


    def _set_distance_from_0(self, dist):
        "changes parameters of plane - the new plane is parallel to old"
        if not self.initialized:
            raise NotInitializedError(None)
        assert dist > 0
        self.p = dist
        self._compute_ABCD_from_angles()


    def get_distance_from_point(self, P):
        "computes distance of plane from point"
        return P * self.cosines - self.p


    def set_as_3points(self, p1, p2, p3):
        "define plane by giving three non-colinear points"
        if len(p1) != 3 or len(p2) != 3 or len(p3) != 3:
            raise UnexpectedArgsError("3 x 3 arguments expected for plane")
        par = lambda r, t: p1[r] * (p2[t] - p3[t]) + p2[r] * (p3[t] - p1[t]) \
                           + p3[r] * (p1[t] - p2[t])
        self.A = par(1,2)
        self.B = par(2,0)
        self.C = par(0,1)
        if (self.A == self.B == self.C == 0):
            raise UnexpectedArgsError("3 colinear points - undefined plane")
        _D = p1[0] * (p2[1]*p3[2] - p3[1]*p2[2]) \
           + p2[0] * (p3[1]*p1[2] - p1[1]*p3[2]) \
           + p3[0] * (p1[1]*p2[2] - p2[1]*p1[2])
        self.D = - _D
        self.initialized = True
        self._compute_angles_from_ABCD()


    def describe_ABCD(self):
        if not self.initialized: return "Plane not initialized"
        return "Plane: %s x + %s y + %s z + %s = 0. " % (self.A, self.B,
                                                         self.C, self.D)

    def describe_angles(self):
        if not self.initialized: return "Plane not initialized"
        return "Plane: x cos(alpha) + y cos(beta) + z cos(gamma) - p = 0. \n" \
                "alpha=%s  beta=%s  gamma=%s  p=%s" % (self.alpha, self.beta,
                                                       self.gamma, self.p)

    def get_normal_vector(self):
        assert self.initialized
        return array((self.A, self.B, self.C), float)

    def get_rotation_matrix_to(self, new_plane):
        n1 = self.get_normal_vector()
        n2 = new_plane.get_normal_vector()
        length_prod = sqrt(inner(n1,n1) * inner(n2,n2))
        cos_angle = dot(n1, n2) / length_prod
        angle = acos(cos_angle)
        p = cross(n1, n2) / length_prod
        sin_angle = sqrt(inner(p,p))
        assert abs(sin_angle - sin(angle)) < 1e-6
        return rodrigues(p / sin_angle, angle, verbose=True)



class LatticePlane(Plane):
    "Plane that can be defined using Miller indices."
    def __init__(self, cell=None, hkl=None, r=None):
        Plane.__init__(self)
        self.set_hkld(hkl, r)
        self.set_cell(cell)

    def __str__(self):
        return self.describe_hkld() + "\n or \n" + Plane.__str__(self)

    def set_hkld(self, hkl, r):
        "define plane using hkl indices and distance from (0,0,0)"
        if hkl is None:
            self.hkl = None
            self.r = r
            return
        if len(hkl) != 3:
            raise UnexpectedArgsError("3 (h k l) arguments for plane"
                                              " required, got %s" % len(hkl))
        self.hkl = array(hkl, float)
        self.r = r
        if not sometrue(self.hkl): #if all are zeros
            raise UnexpectedArgsError("(0 0 0) ???")
        self.compute_plane_parameters()

    def set_cell(self, cell):
        self.cell = cell
        self.compute_plane_parameters()

    def compute_plane_parameters(self):
        if not hasattr(self, "cell") or self.cell is None \
                or not hasattr(self, "hkl") or self.hkl is None:
            return
        #FIXME is there a better way to compute geometrical coordinates
        # of plane? Now this method finds 3 points on plane, ...
        points = []
        for i, index in enumerate(self.hkl):
            if index:
                points.append(self.cell.get_unit_shift(i) / index)
        for i, index in enumerate(self.hkl):
            if not index:
                points.append(points[0] + self.cell.get_unit_shift(i))
        self.set_as_3points(*points)
        self._set_distance_from_0(self.r)


    def describe_hkld(self):
        if not self.initialized: return "Plane not initialized"
        return "Plane indices: (%s %s %s). Distance from origin: %s" % (
                                self.hkl[0], self.hkl[1],self.hkl[2], self.r)


class SurfaceDeformation:
    "Stress-like deformation - perpendicular to surface"
    def __init__(self, depth, fun):
        # depth 0 doesn't prevent deformation, because some atoms can be
        # slightly outside of the surface (i.e. at negative depth)
        self.depth = depth
        self.fun = fun
    def __str__(self):
        return "D_shell=%s" % self.depth


class LatticeSurface(LatticePlane):
    """Surface of CuttedGrain. Stores information about grain shell deformation.
       If hkl is None, it is treated as spherical surface.
    """
    # default SurfaceDeformation
    default_sd = None

    def __init__(self, cell=None, hkl=None, r=None, sd=None):
        LatticePlane.__init__(self, cell, hkl, r)
        self.sd = sd or LatticeSurface.default_sd

    def __str__(self):
        if self.hkl is not None:
            return self.describe_hkld() + ". %s" % self.sd
        elif self.r:
            return "Sphere, r=%s. %s" % (self.r, self.sd)
        else:
            return "Not initialized LatticeSurface"

    def get_planes(self):
        """For sphere - returns bounding planes, eg. cube.
           Otherwise, returns itselfs
        """
        if self.hkl is not None:
            return [self]
        else: #sphere
            r = self.r
            return [Plane((1,0,0,-r)), Plane((0,1,0,-r)), Plane((0,0,1,-r)),
                    Plane((1,0,0,r)), Plane((0,1,0,r)), Plane((0,0,1,r))]



class FreshModel(Model):
    "Base class for configuration generators. Initially Model without atoms."
    def __init__(self, lattice, pbc=None, title=""):
        Model.__init__(self, atoms=[], pbc=pbc, title=title)
        self.lattice = lattice
        self.unit_cell = lattice.unit_cell

    def get_vertices(self):
        assert None, "Virtual method, should be overwritten"

    def compute_scope(self):
        """
        Get minimal and maximal coordinate in a,b,c cell parameters units.
        It cuts out a "rectangular parallelepiped"(?)  that contains the grain.
        """
        self.vertices = self.get_vertices()
        #vertices transformed to orthonormal system
        unit_vertices = [dot(i, self.unit_cell.M) for i in self.vertices]

        normalized = all(node.is_normalized() for node in self.lattice.nodes)
        margin = 0 if normalized else 1

        scopes = []
        for i in zip(*unit_vertices):
            m = floor(round(min(i), 9)) - margin
            M = ceil(round(max(i), 9)) + margin
            scopes.append((int(m), int(M)+1))
        self.a_scope, self.b_scope, self.c_scope = scopes


    def get_scope_info(self):
        "returns info about results of compute_scope() method"
        lar = self.a_scope[1] - self.a_scope[0]
        lbr = self.b_scope[1] - self.b_scope[0]
        lcr = self.c_scope[1] - self.c_scope[0]
        ncl = lar*lbr*lcr
        nnd = len(self.lattice.nodes)
        nat = sum([len(i.atoms_in_node) for i in self.lattice.nodes])
        t = "Considering %ix%ix%i=%i cells, " % (lar, lbr, lcr, ncl)
        t += "%i nodes/cell, %i atoms/cell, " % (nnd, nat)
        t += "%i nodes, %i atoms." % (ncl*nnd, nat*ncl)
        return t


    def get_all_nodes(self):
        "generator of nodes in (and out of) grain"
        for i in range(self.a_scope[0], self.a_scope[1]):
            for j in range(self.b_scope[0], self.b_scope[1]):
                for k in range(self.c_scope[0], self.c_scope[1]):
                    for node in self.lattice.nodes:
                        yield node, array((i, j, k), float) + node.pos_in_cell


    def _do_export_atoms(self, f, format):
        if format == "powdercell":
            self.lattice.export_powdercell(f)
        else:
            Model._do_export_atoms(self, f, format)



class CuttedGrain(FreshModel):
    """Finite grain made from unfinite crystal lattice (CrystalLattice)
       by cleaving using added surfaces (sequence of Plane)"""

    def __init__(self, lattice, surfaces=None, title="generated by gosam"):
        FreshModel.__init__(self, lattice, title=title)
        if isinstance(surfaces, LatticeSurface):
            surfaces = [surfaces]
        for i in surfaces:
            i.set_cell(self.unit_cell)
        self.surfaces = surfaces
        self.operations = []


    def __str__(self):
        return "cutted grain description:\n\n"  + self.lattice.__str__()  \
            + "\n\n%i surfaces:\n  " % len(self.surfaces) \
            + "\n  ".join([i.__str__() for i in self.surfaces]) \
            + "\n\nOperations:\n  " + "\n  ".join(self.operations)


    def export_for_qhull(self):
        planes = []
        for i in self.surfaces:
            planes += i.get_planes()
        s = """\
             3 1
             0      0      0
             4
             %s\n""" % len(planes)
        for i in planes:
            s += "%s %s %s %s\n" % (i.A, i.B, i.C, i.D)
        return s


    def show_in_geomview(self):
        # it shows cubes instead of spheres
        t = self.export_for_qhull()
        p = Popen(["qhull H Fp | qhull G"], shell=True, bufsize=bufsize,
                          stdin=PIPE, stdout=PIPE, close_fds=True)
        #  H - qhalf
        # Fp - print points at halfspace intersections
        #...and compute the convex hull of the intersection points for Geomview
        # G - display Geomview output
        p.stdin.write(t)
        p.stdin.close()
        off_content = p.stdout.read()
        f = NamedTemporaryFile()
        #print "*" * 70
        #print off_content
        f.write(off_content)
        f.flush()
        f.seek(0)
        try:
            commands.getoutput("geomview -nopanels %s" % f.name)
        except KeyboardInterrupt:
            pass


    def get_vertices(self):
        "get vertices of convex hull that contains the grain"
        t = self.export_for_qhull()
        p = Popen(["qhull H Fp"],
                  shell=True, stdin=PIPE, stdout=PIPE)
        p.stdin.write(t)
        p.stdin.close()
        return [[float(j) for j in i.split()]
                    for i in p.stdout.readlines() if len(i.split()) == 3]


    def generate_atoms(self):
        print "generating atoms (atom positions) ..."
        self.atoms = []
        not_included_counter = 0
        shifted_counter = 0
        self.compute_scope()
        print self.get_scope_info()

        #almost inner loop - over all nodes
        for node, abs_pos in self.get_all_nodes():

            node_pos = dot(abs_pos, self.unit_cell.M_1)

            # checking if atoms are inside grain
            outside = False
            for srf in self.surfaces:
                if srf.hkl is not None:
                    t = -(inner(node_pos, srf.cosines) - srf.p)
                else:
                    t = srf.r - sqrt(inner(node_pos, node_pos))
                if t < -1e-12:
                    outside = True
                    not_included_counter += 1
                    break
                srf.tmp_t = t #optimalization ?
            if outside:
                continue

            # inner loop - over (usually one or a few) atoms in node
            for atom in node.atoms_in_node:
                shifts = 0
                min_dist = None
                if atom.non_zero:
                    xyz = dot(abs_pos+atom.pos, self.unit_cell.M_1)
                else:
                    xyz = node_pos
                dxyz = array((0., 0., 0.))
                for srf in self.surfaces:
                    if atom.non_zero:
                        if srf.hkl is not None:
                            t = -(inner(xyz, srf.cosines) - srf.p)
                        else:
                            t = srf.r - sqrt(inner(xyz, xyz))
                    else:
                        t = srf.tmp_t
                    if srf.sd is not None and t < srf.sd.depth:
                        if type(srf.sd.fun) is dict:
                            fun = srf.sd.fun.get(atom.name)
                        else:
                            fun = srf.sd.fun
                        shift = fun(t)
                        if srf.hkl is not None: # plane
                            dxyz -= shift * srf.cosines
                        else: # sphere
                            r0 = srf.r - t
                            dxyz -= shift / r0 * xyz
                        shifts += 1
                    if min_dist is None or t < min_dist:
                        min_dist = t
                self.atoms.append(mdprim.AtomG(atom.name, xyz+dxyz, min_dist))
                if shifts > 0:
                    shifted_counter += 1

        print "%i nodes outside of grain." % not_included_counter
        print "Number of atoms in grain: %i (shifted: %i)" % (len(self.atoms),
                                                              shifted_counter)
        self.log("atom positions generated")




def generate_grain(d):
    """
    Takes dictionary with names from "configuration file", i.e. globals()
    """
    group_nodes =  "do_not_group_nodes" not in d or not d["do_not_group_nodes"]
    if "nodes" in d and "node_atoms" in d:
        if group_nodes:
            nodes = [Node(i, d["node_atoms"]) for i in d["nodes"]]
        else:
            nodes = [Node([atom[1]+node[0], atom[2]+node[1], atom[3]+node[2]],
                          [AtomInNode(atom[0])])
                       for atom in d["node_atoms"] for node in d["nodes"]]
    elif "atoms" in d:
        nodes = [Node(i[1:], [AtomInNode(i[0])]) for i in d["atoms"]]
    else:
        print 'Either "nodes" and "node_atoms" or "atoms" should be defined.'
        return
    lattice = CrystalLattice(d["cell"], nodes)
    g = CuttedGrain(lattice, surfaces=d["surfaces"])
    g.generate_atoms()
    g.print_stochiometry()
    if "vacancy_probability" in d:
        g.make_vacancies(d["vacancy_probability"])
        g.print_stochiometry()
    if "modifier" in d:
        g.modify_atoms(d["modifier"])
    g.round_atom_coordinates()
    if "remove_undercoordinated_atoms" in d:
        bondlength = d["remove_undercoordinated_atoms"]
        g.print_coordination_statistics(bondlength)
        g.remove_undercoordinated_atoms(bondlength)
    extensions = { "pielaszek": "at",
                   "xmol": "xyz",
                   "powdercell": "cel",
                   "dlpoly": "dlpoly",
                   "atomeye": "cfg"
                 }
    def format_name(name):
        if name in extensions:
            return name
        elif name in extensions.values():
            for key, value in extensions.items():
                if value == name:
                    return key
        else:
            msg = "WARNING: unknown file format in 'output_formats': " + name
            logfile.write(msg + "\n")
            print msg

    if "output_formats" in d:
        formats = [format_name(i) for i in d["output_formats"]]
    else:
        formats = ["xmol"]

    # this format requires PBC. We put vacuum (>=1nm) between images
    if "atomeye" in formats:
        g.set_pbc_with_vacuum(width=10)

    if "output_file" in d:
        basename = d["output_file"]
    else:
        filename = os.path.basename(sys.argv[0])
        basename = os.path.splitext(filename)[0]

    for i in formats:
        g.export_atoms(file(basename+"."+extensions[i], 'w'),  format=i)

    logfile = file(basename+".log", 'w')
    logfile.write(str(g) + "\n\n" + "-"*60 + "\n" + file(sys.argv[0]).read())
    return g



if __name__ == '__main__':
    print "Use it as module"

