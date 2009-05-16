# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2
"""\
Atom and CellMethod classes.
"""

from math import sqrt, acos, degrees
import numpy
from numpy import array, inner

import pse
import rotmat

class Atom:
    "atom - symbol and its position in 3D space"
    def __init__(self, name, pos):
        self.name = name
        self.pos = array(pos)

    def __str__(self):
        return "Atom %s at %s" % (self.name, self.pos)

    def get_dist(self, atom2, pbc_half=None):
        d = self.pos - atom2.pos
        if pbc_half is not None and (abs(d) > pbc_half).any():
            for i in range(3):
                if abs(d[i]) > pbc_half[i]:
                    d[i] = 2 * pbc_half[i] - abs(d[i])
        return sqrt(inner(d, d)) # sqrt(sum(d**2)) is slower

    def get_shift(self, atom2, pbc=None):
        d = atom2.pos - self.pos
        if pbc is not None and (abs(d) > pbc / 2).any():
            for i in range(3):
                if d[i] > pbc[i] / 2:
                    d[i] -= pbc[i]
                elif d[i] < -pbc[i] / 2:
                    d[i] += pbc[i]
        return d

    def get_angle(self, atom2, atom3):
        "angle of atom2-self-atom3 in radians"
        v1 = atom2.pos - self.pos
        v2 = atom3.pos - self.pos
        return acos(inner(v1, v2) / sqrt(inner(v1, v1) * inner(v2, v2)))

    def get_temperature(self):
        assert 0, "Unknown temperature"

    def get_velocity(self):
        assert 0, "Unknown velocity"

    def get_ekin(self):
        assert 0, "Unknown velocity"


class AtomG(Atom):
    "atom (in grain) - additional information about distance to surface"
    def __init__(self, name, pos, min_dist):
        Atom.__init__(self, name, pos)
        self.min_dist = min_dist

    def __str__(self):
        return Atom.__str__(self) + "; min. dist: %f" % (self.min_dist)


class AtomVF(Atom):
    "atom with velocity info"
    def __init__(self, name, nr, pos, vel, force):
        Atom.__init__(self, name, pos)
        self.nr = nr
        self.vel = array(vel)
        self.force = array(force)

    def __str__(self):
        return "Atom %s at %s; velocity: %s" % (self.name, self.pos, self.vel)

    def get_mass(self):
        return pse.get_atom_mass(self.name)

    def get_velocity(self):
        return sqrt(inner(self.vel, self.vel))

    def get_ekin(self):
        conv_factor = 1.0364269e-10 # u * (A/ns)^2 -> eV
        return conv_factor * 0.5 * self.get_mass() * inner(self.vel, self.vel)

    def get_temperature(self):
        if self.vel is None or self.vel[0] is None:
            return 0
        kB = 8.617343e-5 # eV/K
        #kB = 831447.2 # u A^2 / ns^2 K
        return 2 * self.get_ekin() / (3*kB)


#only supports orthorhombic PBC
class CellMethod:
    "Neighbour list maker"
    def __init__(self, atoms, r, pbc=None):
        self.atoms = atoms
        if pbc is not None:
            assert rotmat.is_diagonal(pbc)
            self.box_d = self.pbc = pbc.diagonal()
        else:
            box_min, box_max = self._find_containing_box()
            self.box_d = box_max - box_min
            self.pbc = None
        self._make_cells(r)

    def _find_containing_box(self):
        m = M = self.atoms[0].pos
        for i in self.atoms:
            m = numpy.minimum(i.pos, m)
            M = numpy.maximum(i.pos, M)
        return m, M

    def _get_cell_coord(self, a):
        return numpy.floor(a.pos * self.nc / self.box_d).astype(int) % self.nc

    def _make_cells(self, r):
        self.r = r
        self.nc = (self.box_d / (r + 1e-9)).astype(int)

        # avoid memory problems when cut-off is small
        nc_prod = self.nc[0] * self.nc[1] * self.nc[2]
        if nc_prod > 1000000 and nc_prod > 2 * len(self.atoms):
            f = float(nc_prod) / len(self.atoms)
            self.nc /= f**(1./3)

        cell_count = self.nc[0] * self.nc[1] * self.nc[2]
        self.cells = [[] for i in range(cell_count)]
        for a_idx, a in enumerate(self.atoms):
            nx, ny, nz = self._get_cell_coord(a)
            cell_idx = (nx * self.nc[1] + ny) * self.nc[2] + nz
            self.cells[cell_idx].append(a_idx)
        print "... system divided into %d x %d x %d cells ..." % tuple(self.nc)


    def _get_neigh_cells_in_dim(self, n, c):
        assert 0 <= n < c
        if c > 2:
            return (n-1) % c, n, (n+1) % c
        elif c == 2:
            return 0, 1
        elif c == 1:
            return 0,

    def _get_neighbour_cells(self, a):
        nx, ny, nz = self._get_cell_coord(a)
        for tix in self._get_neigh_cells_in_dim(nx, self.nc[0]):
            for tiy in self._get_neigh_cells_in_dim(ny, self.nc[1]):
                for tiz in self._get_neigh_cells_in_dim(nz, self.nc[2]):
                    cell_idx = (tix * self.nc[1] + tiy) * self.nc[2] + tiz
                    yield cell_idx


    def get_neighbours(self, n, extra_condition=None):
        assert self.pbc is None
        a = self.atoms[n]
        for cell_idx in self._get_neighbour_cells(a):
            for i in self.cells[cell_idx]:
                dist = a.get_dist(self.atoms[i])
                if n != i and dist < self.r and (extra_condition is None
                                                    or extra_condition(dist)):
                    yield i


    def pop_neighbours(self, n, max_dist=None):
        if max_dist is None:
            max_dist = self.r
        if self.pbc is not None:
            half_pbc = array(self.pbc) / 2.
        else:
            half_pbc = None
        a = self.atoms[n]
        for cell_idx in self._get_neighbour_cells(a):
            cell = self.cells[cell_idx]
            for i_idx, i in enumerate(cell):
                if n != i and (a.get_dist(self.atoms[i], half_pbc) < max_dist):
                    yield cell.pop(i_idx)


    def count_neighbours(self, n, extra_condition=None):
        counter = 0
        for i in self.get_neighbours(n, extra_condition):
            counter += 1
        return counter


    def get_atoms_to_remove(self):
        """Return map, which keys are indices of atoms for removing,
           and values are neighbours that caused atom to be removed.
           Atoms are removed to leave not more than one atom in one cell.
        """
        to_be_deleted = {}
        for n, i in enumerate(self.atoms):
            for j in self.pop_neighbours(n):
                if j in to_be_deleted:
                    if n not in to_be_deleted[j]:
                        to_be_deleted[j].append(n)
                elif n in to_be_deleted:
                    if j not in to_be_deleted[n]:
                        to_be_deleted[n].append(j)
                else:
                    to_be_deleted[j] = [n]
        return to_be_deleted


