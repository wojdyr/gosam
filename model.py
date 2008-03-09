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

class Model:
    """Configuration -- atoms, PBC (only parallelepiped)"""
    def __init__(self, atoms, pbc, title=""):
        self.atoms = atoms
        self.pbc = pbc 
        self.title = title

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


    def _get_atoms_to_remove(self, cm):
        """return map, which keys are indices of atoms for removing,
           and values are neighbours that caused atom to be removed
        """
        to_be_deleted = {}
        for n, i in enumerate(self.atoms):
            for j in cm.pop_neighbours(n):
                if j in to_be_deleted:
                    if n not in to_be_deleted[j]:
                        to_be_deleted[j].append(n)
                elif n in to_be_deleted:
                    if j not in to_be_deleted[n]:
                        to_be_deleted[n].append(j)
                else:
                    to_be_deleted[j] = [n]
        return to_be_deleted
                

    def _print_deleted_dist_stats(self, to_be_deleted):
        dd = []
        pbc_half = array(self.pbc.diagonal()) / 2.
        for k,v in to_be_deleted.iteritems():
            for j in v:
                dist = self.atoms[k].get_dist(self.atoms[j], pbc_half=pbc_half)
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


    def remove_close_neighbours(self, distance):
        """Remove atoms in such a way that no two atoms are in distance
        smaller than `distance'
        """
        assert rotmat.is_diagonal(self.pbc)
        print "Removing atoms in distance < %f ..." % distance
        before = len(self.atoms)
        cm = mdprim.CellMethod(self.atoms, distance, self.pbc)
        to_be_deleted = self._get_atoms_to_remove(cm)
        self._print_deleted_dist_stats(to_be_deleted)
        #self._shift_before_removing(to_be_deleted)
        tbd_idx = to_be_deleted.keys()
        tbd_idx.sort(reverse=True)
        for i in tbd_idx:
            del self.atoms[i]
        rem = before - len(self.atoms)
        print "... %i atoms removed." % rem
        self.print_stochiometry()
        self.log("removed %i too-close atoms" % rem)


    def export_atoms(self, f, format=None):
        """
        save atoms to file f in one of possible formats 
        """
        if type(f) in (str, unicode):
            f = file(f, 'w')
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



