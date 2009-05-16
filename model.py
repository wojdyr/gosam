# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2
"""\
class Model -- atomistic model, optionally in PBC -- vacancies,
colision detection, saving to file, etc.
"""

import inspect
import random
import sys
import numpy
from numpy import array

import mdprim
import mdfile
import rotmat


def _sort_and_uniq(dd):
    "sort real number and leave only unique ones (compare using epsilon)"
    dd.sort()
    n = 0
    while n+1 < len(dd):
        if abs(dd[n] - dd[n+1]) < 1e-6:
            del dd[n+1]
        else:
            n += 1


class Model:
    """Configuration -- atoms, PBC (only parallelepiped)"""
    def __init__(self, atoms, pbc, title="", comments=None):
        self.atoms = atoms
        self.pbc = pbc
        self.title = title
        self.comments = comments

    def log(self, text):
        "to enable logging, do 'operations=[]'"
        if hasattr(self, "operations"):
            self.operations.append(text)


    def round_atom_coordinates(self, ndigits=9):
        "rounds coordinates of all atoms (points)"
        for i in self.atoms:
            for j in range(3):
                i.pos[j] = round(i.pos[j], ndigits) + 0.0
        self.log("atom coordinates rounded to %s digits after dot." % ndigits)


    def make_vacancies(self, vacancy_probability):
        " vacancy_probability: {atom_name: float <0,1>} "
        if not vacancy_probability:
            return
        len_before = len(self.atoms)
        for idx, atom in enumerate(self.atoms):
            if type(vacancy_probability) is dict:
                try:
                    p = vacancy_probability[atom.name]
                except KeyError:
                    continue
            else: # vacancy_probability is function
                p = vacancy_probability(atom)
            if random.random() < p:
                #print idx, atom
                del self.atoms[idx]
        len_after = len(self.atoms)
        n_vacancies = len_before-len_after
        print "%i atoms were deleted (vacancies were made). %i atoms left." % (
                                                        n_vacancies, len_after)
        self.log("vacancies where generated (%i of %i) using probabilities: %s"
                   % (n_vacancies, len_before, vacancy_probability))


    def modify_atoms(self, modifier):
        "Use given function to modify atom (point) coordinates"
        if not modifier:
            return
        for i in self.atoms:
            modifier(i)
        self.log("atom positions modified using function: \n\t"
                           + inspect.getsource(modifier).replace("\n", "\n\t"))


    def count_neighbours(self, atom, max_bondlength):
        "O(N^2); use mdprim.CellMethod instead"
        print "WARNING: ineffective neighbour counting in use"
        neighbors = 0
        for j in self.atoms:
            #optimization
            if abs(j.pos[0] - atom.pos[0]) < max_bondlength \
                    and abs(j.pos[1] - atom.pos[1]) < max_bondlength \
                    and 1e-3 < atom.get_dist(j) < max_bondlength:
                neighbors += 1
        return neighbors


    def print_coordination_statistics(self, max_bondlength):
        print "Coordination statistics: ",
        sys.stdout.flush()
        stat = {}
        #for i in self.atoms:
        #    n = self.count_neighbours(i, max_bondlength)
        #    stat[n] = stat.get(n, 0) + 1
        cm = mdprim.CellMethod(self.atoms, max_bondlength)
        for a_idx, i in enumerate(self.atoms):
            n = cm.count_neighbours(a_idx)
            stat[n] = stat.get(n, 0) + 1
        s = stat.items()
        s.sort(lambda x,y: -cmp(x[1], y[1]))
        print ", ".join("%i: %i" % i for i in s)


    def print_stochiometry(self):
        print "Stochiometry: ",
        stat = {}
        for i in self.atoms:
            if i.name in stat:
                stat[i.name] += 1
            else:
                stat[i.name] = 1
        s = stat.items()
        s.sort(lambda x,y: -cmp(x[1], y[1]))
        print ", ".join("%s: %i" % i for i in s)


    def remove_undercoordinated_atoms(self, max_bondlength):
        """Remove atoms that have only 1 nearest neighbor
        and some with 2 nearest neighbors (stoichiometry is conserved).
        For use in tetrahedrally coordinated lattices.
        """
        print "Removing under-coordinated atoms..."
        before = len(self.atoms)
        for iter in range(2):
            cm = mdprim.CellMethod(self.atoms, max_bondlength)
            to_be_deleted = []
            for n, i in enumerate(self.atoms):
                c = cm.count_neighbours(n)
                #c = self.count_neighbours(i, max_bondlength)
                if c <= 1:
                    to_be_deleted.append(n)
            to_be_deleted.sort(reverse=True)
            for i in to_be_deleted:
                del self.atoms[i]
        rem = before - len(self.atoms)
        print "... %i atoms removed." % rem
        self.print_coordination_statistics(max_bondlength)
        self.print_stochiometry()
        self.log("removed %i under-coordinated atoms" % rem)


    def _print_deleted_dist_stats(self, atoms, to_be_deleted):
        dd = []
        pbc_half = array(self.pbc.diagonal()) / 2.
        for k,v in to_be_deleted.iteritems():
            for j in v:
                dist = atoms[k].get_dist(atoms[j], pbc_half=pbc_half)
                dd.append(dist)
        if not dd:
            print "no atoms were too close"
            return
        print "   deleted atoms distances: from %s to %s" % (min(dd), max(dd))


    def _shift_before_removing(self, to_be_deleted):
        """if only pairs of atoms of the same species are too close
           to each other, move the atom that won't be deleted to the position
           between it's old position and the position of the neighbour.
        """
        pbc = self.pbc.diagonal()
        for k, v in to_be_deleted.iteritems():
            assert len(v) == 1, "%s %s" % (k, v)
            assert v[0] not in to_be_deleted, "%s" % v[0]
            a = self.atoms[k]
            b = self.atoms[v[0]]
            assert a.name == b.name, "%s %s" % (a.name, b.name)
            # a will be deleted, b not
            d = b.get_shift(a, pbc=pbc)
            b.pos += d / 2


    def get_atoms_to_be_removed(self, atoms, distance):
        assert rotmat.is_diagonal(self.pbc)
        cm = mdprim.CellMethod(atoms, distance, self.pbc)
        return cm.get_atoms_to_remove()

    def remove_close_neighbours(self, distance, atoms=None):
        """Remove atoms in such a way that no two atoms are in distance
        smaller than `distance'
        """
        if atoms is None:
            atoms = self.atoms
        before = len(atoms)
        to_be_deleted = self.get_atoms_to_be_removed(atoms, distance)
        self._print_deleted_dist_stats(atoms, to_be_deleted)
        #self._shift_before_removing(to_be_deleted)
        tbd_idx = to_be_deleted.keys()
        tbd_idx.sort(reverse=True)
        for i in tbd_idx:
            del atoms[i]
        rem = before - len(atoms)
        print "... %i atoms removed." % rem
        if atoms is self.atoms: # otherwise self.atoms stats are useless
            self.print_stochiometry()
        self.log("removed %i too-close atoms" % rem)


    def add_close_neigh_properties(self):
        """
        add .r1 and .r2 members (None or float) to each atom
        r1 - the closest distance to other atom
        r2 - the closest distance to other atom with the same symbol
        """
        r1_max = 1.87
        r2_max = 3.00

        for i in self.atoms:
            i.r1 = None
            i.r2 = None

        pbc_half = array(self.pbc.diagonal()) / 2.

        # r1
        to_be_rm1 = self.get_atoms_to_be_removed(self.atoms, r1_max)
        for k,v in to_be_rm1.iteritems():
            atom = self.atoms[k]
            d = min(atom.get_dist(self.atoms[j], pbc_half=pbc_half) for j in v)
            atom.r1 = d

        # r2
        a_name = self.atoms[0].name
        a_atoms = [i for i in self.atoms if i.name == a_name]
        b_atoms = [i for i in self.atoms if i.name != a_name]
        for x_atoms in a_atoms, b_atoms:
            to_be_rm2 = self.get_atoms_to_be_removed(x_atoms, r2_max)
            for k,v in to_be_rm2.iteritems():
                atom = x_atoms[k]
                d = min(atom.get_dist(x_atoms[j], pbc_half=pbc_half) for j in v)
                #if not atom.r1 or d > atom.r1:
                atom.r2 = d


    def output_all_removal_possibilities(self, filename):
        assert "%" in filename

        for n, i in enumerate(self.atoms):
            i.nr = n

        self.add_close_neigh_properties()

        distances1 = [0] + [i.r1 + 1e-6 for i in self.atoms if i.r1]
        distances2 = [0] + [i.r2 + 1e-6 for i in self.atoms if i.r2]

        _sort_and_uniq(distances1)
        _sort_and_uniq(distances2)
        print "inter-atomic distances:", distances1
        print "same species distances:", distances2
        print "atoms count:", len(self.atoms)
        print

        counter = 1
        orig_atoms = self.atoms
        all_rm = []
        for j in distances2:
            for i in distances1:
                #if j <= i:
                #    continue
                rm = [a.nr for a in orig_atoms if (a.r1 and a.r1 < i)
                                                   or (a.r2 and a.r2 < j)]
                if rm in all_rm:
                    print "ignore cutoffs: %g, %g (%d atoms)" % (i, j, len(rm))
                    continue
                all_rm.append(rm)
                self.title = "del %d atoms with cutoffs: %g, %g" % (
                                                                len(rm), i, j)
                print self.title
                self.atoms = [a for a in orig_atoms if a.nr not in rm]
                fn = filename.replace('%', str(counter))
                self.export_atoms(fn)
                counter += 1


    def _find_symmetric_z_distances(self):
        distances = [0]
        for i in self.atoms:
            d = 2 * i.pos[2]
            if 1e-7 < d < 3.0:
                distances.append(d + 1e-6)
        _sort_and_uniq(distances)
        print "same species distances:", distances
        print "atoms count:", len(self.atoms)
        return distances

    def output_all_removal2_possibilities(self, filename):
        assert "%" in filename
        distances = self._find_symmetric_z_distances()
        orig_atoms = self.atoms
        for n, j in enumerate(distances):
            self.atoms = [a for a in orig_atoms if not 1e-7 < a.pos[2] < j / 2.]
            ndel = len(orig_atoms) - len(self.atoms)
            self.title = "del %d atoms with cutoff: %g" % (ndel, j)
            print self.title
            fn = filename.replace('%', str(n))
            self.export_atoms(fn)


    def output_all_possibilities_all_stoich(self, filename):
        assert "%" in filename
        distances = self._find_symmetric_z_distances()
        orig_atoms = self.atoms
        species = self.count_species()
        assert len(species) == 2
        name1, name2 = sorted(species.keys()) # C, Si
        for n1, j1 in enumerate(distances):
            for n2, j2 in enumerate(distances):
                zmax = { name1: j1 / 2., name2: j2 / 2. }
                self.atoms = [a for a in orig_atoms
                              if not 1e-7 < a.pos[2] < zmax[a.name]]
                ndel = len(orig_atoms) - len(self.atoms)
                self.title = "del %d atoms with cutoffs: %g, %g" % (ndel,j1,j2)
                print self.title
                fn = filename.replace('%', "%d-%d" % (n1, n2))
                self.export_atoms(fn)


    def export_atoms(self, f, format=None):
        """
        save atoms to file f in one of possible formats
        """
        if type(f) in (str, unicode):
            f = mdfile.open_any(f, 'w')
        if format is None:
            format = mdfile.get_type_from_filename(f.name);
        format = format.lower()
        print "Saving atoms to file '%s' in format '%s'" % (f.name, format)
        self._do_export_atoms(f, format)
        self.log("atoms saved to file '%s' in format '%s'" % (f.name, format))

    def _do_export_atoms(self, f, format):
        if format == "xmol":
            mdfile.export_as_xmol(self.atoms, f, self.title)
        elif format == "pielaszek":
            mdfile.export_for_pielaszek(self.atoms, f)
        elif format == "dlpoly":
            mdfile.export_for_dlpoly(self.atoms, f, self.title)
        elif format == "atomeye":
            mdfile.export_for_atomeye(self, f)
        elif format == "poscar":
            mdfile.export_as_poscar(self, f)
        elif format == "gulp":
            mdfile.export_as_gulp(self, f)
        elif format == "lammps":
            mdfile.export_as_lammps(self, f)
        else:
            print >>f, "Unknown format requested: %s" % format


    def get_center(self, onAtom=False):
        n = len(self.atoms)
        ctr_pos = sum([i.pos for i in self.atoms]) / n
        ctr = mdprim.Atom("<Center>", ctr_pos)
        if onAtom: # the nearest atom
            dists = [ctr.get_dist(i) for i in self.atoms]
            return self.atoms[dists.index(min(dists))]
        else: # the center
            return ctr


    def get_T_vs_centerdist(self, n=100):
        ctr = self.get_center()
        t = [(ctr.get_dist(i), i.get_temperature()) for i in self.atoms]
        print "Average temperature:", sum(i[1] for i in t) / len(t),
        print "max.", max(i[1] for i in t)
        t.sort(lambda x,y: cmp(x[0], y[0]))
        xy = []
        for i in range(len(t)//n):
            g = t[i*n: (i+1)*n]
            x = sum(j[0] for j in g) / n
            y = sum(j[1] for j in g) / n
            xy.append((x, y))
        return xy


    def write_T_vs_centerdist(self, filename, n_group=100):
        print "Writing radial distribution of temperature " \
                "(atoms grouped %s) to file %s" % (n_group, filename)
        ofile = file(filename, 'w')
        for i in self.get_T_vs_centerdist(n_group):
            print >>ofile, i[0], i[1]

    def set_pbc_with_vacuum(self, width):
        pbc = numpy.zeros((3,3))
        for i in range(3):
            k = lambda atom: atom.pos[i]
            pbc[i][i] = k(max(self.atoms, key=k)) - k(min(self.atoms, key=k)) \
                        + width
        self.pbc = pbc


    def count_species(self):
        counts = {}
        for i in self.atoms:
            if i.name in counts:
                counts[i.name] += 1
            else:
                counts[i.name] = 1
        return counts


