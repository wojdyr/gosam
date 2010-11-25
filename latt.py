# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2
"""\
CrystalLattice, UnitCell, etc.
"""

from math import cos, sin, acos, asin, sqrt, pi, floor, ceil, radians, degrees
from numpy import linalg
from numpy import array, dot, transpose


class UnitCell:
    """basic unit cell  - triclinic
    """
    def __init__(self, a, b, c, alpha, beta, gamma,
                 system="triclinic", rad=False, reciprocal=None):
        self.abc = self.a, self.b, self.c = a, b, c
        if not rad: #degrees
            self.alpha, self.beta, self.gamma = alpha, beta, gamma
            self.alpha_rad, self.beta_rad, self.gamma_rad = radians(alpha), \
                                                radians(beta), radians(gamma)
        else: #radians
            self.alpha, self.beta, self.gamma = degrees(alpha), \
                                                degrees(beta), degrees(gamma)
            self.alpha_rad, self.beta_rad, self.gamma_rad = alpha, beta, gamma
        self._compute_sin_cos_V()
        if not reciprocal: #lattice in direct space
            self.is_reciprocal = False
            self.reciprocal = self.get_reciprocal_unit_cell()
        else: #lattice in reciprocal space
            self.is_reciprocal = not reciprocal.is_reciprocal
            self.reciprocal = reciprocal
        self.system = system
        self.system_name = self.system.title() + " system " \
                         + (self.is_reciprocal and "(reciprocal)" or "(direct)")
        self.compute_transformation_matrix()


    def __str__(self):
        return self.system_name \
                + " a=%s, b=%s, c=%s, alpha=%s, beta=%s, gamma=%s" % (self.a,
                            self.b, self.c, self.alpha, self.beta, self.gamma)

    # overloaded for hexagonal lattice
    def get_orthorhombic_supercell(self):
        if self.alpha == 90 and self.beta == 90 and self.gamma == 90:
            return self.a, self.b, self.c
        else:
            assert "Not implemented."

    def _compute_sin_cos_V(self):
        "precomputations of some values (eg. sin(alpha)) -- for optimization"
        a, b, c = self.a, self.b, self.c
        alpha, beta, gamma = self.alpha_rad, self.beta_rad, self.gamma_rad
        self.sines = array([sin(alpha), sin(beta), sin(gamma)])
        self.cosines = array([cos(alpha), cos(beta), cos(gamma)])
        #Giacovazzo p.62
        G = (a*b*c)**2 * (1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2
                          + 2*cos(alpha)*cos(beta)*cos(gamma))
        self.V = sqrt(G)


    def get_reciprocal_unit_cell(self):
        """ returns instance of UnitCell, that is reciprocal to self
        self.reciprocal.(a|b|c|alpha|...)
            a* = self.reciprocal.a, etc.
            (Giacovazzo, p. 64)
        """
        a, b, c, V = self.a, self.b, self.c, self.V
        alpha, beta, gamma = self.alpha_rad, self.beta_rad, self.gamma_rad
        ar = b*c*sin(alpha)/V
        br = a*c*sin(beta)/V
        cr = a*b*sin(gamma)/V
        cos_alphar = (cos(beta)*cos(gamma)-cos(alpha)) / (sin(beta)*sin(gamma))
        cos_betar = (cos(alpha)*cos(gamma)-cos(beta)) / (sin(alpha)*sin(gamma))
        cos_gammar = (cos(alpha)*cos(beta)-cos(gamma)) / (sin(alpha)*sin(beta))
        return UnitCell(ar, br, cr, acos(cos_alphar), acos(cos_betar),
                              acos(cos_gammar), rad=True, reciprocal=self)


    def compute_transformation_matrix(self):
        """ sets self.M and self.M_1 (M^-1)
             [       1/a               0        0   ]
         M = [  -cosG/(a sinG)    1/(b sinG)    0   ]
             [     a*cosB*         b*cosA*      c*  ]

         where A is alpha, B is beta, G is gamma, * means reciprocal space.
         E=MA
         (Giacovazzo, p.68)
        """
        a, b, c = self.a, self.b, self.c
        alpha, beta, gamma = self.alpha_rad, self.beta_rad, self.gamma_rad
        r = self.reciprocal
        self.M = array([
             (1/a,  0,  0),
             (-cos(gamma)/(a * sin(gamma)),  1/(b*sin(gamma)),  0),
             (r.a * cos(r.beta_rad),  r.b * cos(r.alpha_rad),  r.c)
        ])

        self.M_1 = array([
             (a,  0,  0),
             (b*cos(gamma),  b*sin(gamma),  0),
             (c*cos(beta),  -c*sin(beta)*cos(r.alpha_rad),  1/r.c)
        ])

    def rotate(self, rot_mat):
        assert (abs(transpose(rot_mat) - linalg.inv(rot_mat)) < 1e-6).all(), \
                "not orthogonal"
        assert abs(linalg.det(rot_mat) - 1) < 1e-6, "not a pure rotation matrix"
        self.M_1 = dot(self.M_1, rot_mat)
        self.M = dot(transpose(rot_mat), self.M)
        #print "1=", dot(self.M_1, self.M)

    def get_unit_shift(self, i):
        return self.M_1[i]



class CubicUnitCell(UnitCell):
    def __init__(self, a):
        UnitCell.__init__(self, a, a, a, 90, 90, 90, system="cubic")

    def __str__(self):
        return self.system_name + "  a=%s" % self.a


class TetragonalUnitCell(UnitCell):
    def __init__(self, a, c):
        UnitCell.__init__(self, a, a, c, 90, 90, 90, system="tetragonal")

    def __str__(self):
        return self.system_name + "  a=%s, c=%s" % (self.a, self.c)


class OrthorhombicUnitCell(UnitCell):
    def __init__(self, a, b, c):
        UnitCell.__init__(self, a, b, c, 90, 90, 90, system="orthorhombic")

    def __str__(self):
        return self.system_name + "  a=%s, b=%s, c=%s"%(self.a, self.b, self.c)


class HexagonalUnitCell(UnitCell):
    def __init__(self, a, c):
        UnitCell.__init__(self, a, a, c, 90, 90, 120, system="hexagonal")

    def __str__(self):
        return self.system_name + "  a=%s, c=%s" % (self.a, self.c)

    def get_orthorhombic_supercell(self):
        return self.a, self.a * 3**0.5, self.c



class AtomInNode:
    """atom with its coordinates (as fraction of unit cell parameters) in node
      (atoms in one node can't be separated; node is in unit cell) """
    def __init__(self, name, xa=0, yb=0, zc=0):
        self.name = name
        self.pos = array((xa, yb, zc), float)
        self.non_zero = (xa != 0 or yb != 0 or zc != 0)

    def __str__(self):
        return "%s at %s in node"% (self.name, tuple(self.pos))


class Node:
    """ Node in unit cell consists of a few (eg. 1 or 2) atoms, which should
        be kept together when cutting grains
    """
    def __init__(self, pos_in_cell, atoms_in_node):
        self.pos_in_cell = array(pos_in_cell, float)
        self.atoms_in_node = [isinstance(i, AtomInNode) and i or AtomInNode(*i)
                              for i in atoms_in_node]

    def __str__(self):
        return "Node at %s in cell with: %s" % (tuple(self.pos_in_cell),
                               ", ".join([str(i) for i in self.atoms_in_node]))

    def shift(self, v):
        self.pos_in_cell = (self.pos_in_cell + array(v)) % 1.0

    def is_normalized(self):
        """Are positions of all the atoms in this node in the <0,1) range.
           (i.e. checking self.pos_in_cell+atom.pos)
        """
        for atom in self.atoms_in_node:
            p = self.pos_in_cell + atom.pos
            if (p >= 1).any() or (p < 0).any():
                return False
        return True


class CrystalLattice:
    "Unfinite crystal lattice = unit cell + nodes in cell"
    def __init__(self, unit_cell, nodes, name="Mol"):
        self.unit_cell = unit_cell
        self.nodes = nodes
        assert isinstance(name, basestring)
        self.name = name

    def __str__(self):
        return str(self.unit_cell) + "   Nodes:\n" \
                + "\n".join([str(i) for i in self.nodes])

    def count_species(self):
        names = set()
        for i in self.nodes:
            for j in i.atoms_in_node:
                names.add(j.name)
        return len(names)

    def swap_node_atoms_names(self):
        "works only if there are two atoms in every node"
        for i in self.nodes:
            ain = i.atoms_in_node
            assert len(ain) == 2
            ain[0].name, ain[1].name = ain[1].name, ain[0].name

    def shift_nodes(self, v):
        for i in self.nodes:
            i.shift(v)

    def export_powdercell(self, f):
        cell = self.unit_cell
        print >>f, "CELL %f %f %f %f %f %f" % (cell.a, cell.b, cell.c,
                                             cell.alpha, cell.beta, cell.gamma)
        names = []
        for node in self.nodes:
            for atom in node.atoms_in_node:
                pos = node.pos_in_cell + atom.pos
                names.append(atom.name)
                xname = "%s%i" % (atom.name, names.count(atom.name))
                print >>f, "%-4s %-4s %-8s %-8s %-8s" % (atom.name, atom.name,
                                                        pos[0], pos[1], pos[2])


def generate_polytype(a, h, polytype):
    """utility for dealing with various polytypes. a is a from hexagonal
       structure, h is a distance between layers.
       Usage:
         cell, nodes = generate_polytype(a=3.073, h=2.51, polytype="ABABABC")
    """
    pos = { "A" : (0.0, 0.0),
            "B" : (1/3., 2/3.),
            "C" : (2/3., 1/3.)
          }
    polytype = polytype.upper()
    N = len(polytype)
    cell = HexagonalUnitCell(a=a, c=h*N)
    nodes = []
    for n, i in enumerate(polytype):
        t = pos[i]
        nodes.append((t[0], t[1], float(n)/N))
    return cell, nodes



