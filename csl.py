#!/usr/bin/env python
# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2
"""\
Coincidence Site Lattice related utilities.
"""
usage_string = """\
 Usage:
  csl.py hkl [limit=M] [max_angle=A] - list all CSL boundaries with rotation
                                       axis [hkl], sigma < M (default: 1000)
                                       and angle < A (default: 90)

  csl.py hkl sigma     - details about CSL with given sigma

  csl.py hkl m n       - details about CSL generated with integers m and n

 Examples:
  csl.py 100                            # generates a very long list of GBs
  csl.py 100 limit=50 max_angle=45      # generates a shorter list
  csl.py 111 31                         # show details of the specified GB
"""

import sys
import functools
from math import degrees, atan, sqrt, pi, ceil
import numpy
from numpy import array, identity, dot, inner
from numpy.linalg import inv, det, solve
from rotmat import rodrigues, print_matrix

def gcd(a, b):
    "Returns the Greatest Common Divisor"
    assert isinstance(a, int)
    assert isinstance(b, int)
    while a:
        a, b = b%a, a
    return b

# 0 is coprime only with 1
def coprime(a, b):
    return gcd(a,b) in (0, 1)

def gcd_array(a):
    r = abs(a[0])
    for i in a[1:]:
        if i != 0:
            r = gcd(r, abs(i))
    return r


def parse_miller(s):
    if len(s) == 3 and s.isdigit():
        return array([int(s[i]) for i in range(3)])
    elif ',' in s:
        sp = s.split(",")
        assert len(sp) == 3
        return array([int(i) for i in sp])
    else:
        raise ValueError("Can't parse miller indices: %s" % s)


def get_cubic_sigma(hkl, m, n=1):
    sqsum = inner(hkl, hkl)
    sigma = m*m + n*n * sqsum
    while sigma != 0 and sigma % 2 == 0:
        sigma /= 2
    return (sigma if sigma > 1 else None)

def get_cubic_theta(hkl, m, n=1):
    h,k,l = hkl
    sqsum = h*h + k*k + l*l
    assert sqsum > 0
    if m > 0:
        return 2 * atan(sqrt(sqsum) * n / m)
    else:
        return pi

def get_theta_m_n_list(hkl, sigma, verbose=False):
    if sigma == 1:
        return [(0., 0, 0)]
    thetas = []

    # From Grimmer, Acta Cryst. (1984). A40, 108-112
    #    S = m^2 + (u^2+v^2+w^2) n^2     (eq. 2)
    #    S = alpha * Sigma               (eq. 4)
    #   where alpha = 1, 2 or 4.
    # Since (u^2+v^2+w^2) n^2 > 0,
    # thus alpha * Sigma > m^2    =>   m^2 < 4 * Sigma
    max_m = int(ceil(sqrt(4*sigma)))

    for m in range(max_m):
        for n in range(1, max_m):
            if not coprime(m, n):
                continue
            s = get_cubic_sigma(hkl, m, n)
            if s != sigma:
                continue
            theta = get_cubic_theta(hkl, m, n)
            if verbose:
                print "m=%i n=%i" % (m, n), "%.3f" % degrees(theta)
            thetas.append((theta, m, n))
    return thetas

def find_theta(hkl, sigma, verbose=True, min_angle=None):
    thetas = get_theta_m_n_list(hkl, sigma, verbose=verbose)
    if min_angle:
        thetas = [i for i in thetas if i[0] >= min_angle]
    if thetas:
        return min(thetas, key= lambda x: x[0])


def _get_unimodular_transformation():
    "generator of a few possible unimodular transformations"
    # randomly choosen
    yield identity(3)
    yield array([[1, 0, 1],
                [0, 1, 0],
                [0, 1, 1]])
    yield array([[1, 0, 1],
                [0, 1, 0],
                [0, 1, -1]])
    yield array([[1, 0, 1],
                [0, 1, 0],
                [-1,1, 0]])
    yield array([[1, 0, 1],
                 [1, 1, 0],
                 [1, 1, 1]])


def _get_S():
    Sp = identity(3) # for primitive cubic
    Sb = array([[0.5, -0.5, 0],
                [0.5,  0.5, 0],
                [0.5,  0.5, 1]]) # body-centered cubic
    Sf = array([[0.5, 0.5, 0],
                [0,   0.5, 0.5],
                [0.5, 0,   0.5]]) # face-centered cubic
    # Sf doesn't work?
    return Sp


def transpose_3x3(f):
    """decorator; transpose the first argument and the return value (both
    should be 3x3 arrays). This makes column operations easier"""
    @functools.wraps(f)
    def wrapper(*args, **kwds):
        args_list = list(args)
        assert args_list[0].shape == (3,3)
        args_list[0] = args_list[0].transpose()
        ret_val = f(*args_list, **kwds)
        assert ret_val.shape == (3,3)
        return ret_val.transpose()
    return wrapper


@transpose_3x3
def beautify_matrix(T):
    # We don't want to change the lattice.
    # We use only elementary column operations that don't change det
    def looks_better(a, b):
        x = numpy.abs(a)
        y = numpy.abs(b)
        #return x.sum() < y.sum()
        #return x.sum() < y.sum() or (x.sum() == y.sum() and x.max() < y.max())
        return x.max() < y.max() or (x.max() == y.max() and x.sum() < y.sum())

    def try_add(a, b):
        changed = False
        while looks_better(a+b, a):
            a += b
            changed = True
        return changed

    def try_add_sub(a, b):
        return try_add(a, b) or try_add(a, -b)

    while True:
        changed = False
        for i in range(3):
            for j in range(3):
                if i != j and not changed:
                    changed = try_add_sub(T[i], T[j])
                    if changed:
                        break
        if not changed:
            break

    return T


@transpose_3x3
def make_parallel_to_axis(T, col, axis):
    """\
    T: matrix 3x3, i.e. 3 vectors, 2*T is integer matrix
    axis: vector (3)
    return value:
       matrix T is transformed using operations:
         - interchanging two columns
         - adding a multiple of one column to another,
         - multiplying column by -1
       such that the result matrix has the same det
                                and has first vector == axis
       the transformation is _not_ rotation
    """
    double_T = False
    if not is_integer(T):
        T *= 2 # make it integer, will be /=2 at the end
        double_T = True
    axis = array(axis)
    c = solve(T.transpose(), axis) # c . T == axis
    if not is_integer(c):
        mult = find_smallest_multiplier(c)
        c *= mult
    c = c.round().astype(int)
    #print "c", c
    sel_val = min([i for i in c if i != 0], key=abs)
    if abs(sel_val) != 1: # det must be changed
        print "\n\tWARNING: Volume increased by %i" % abs(sel_val)
    idx = c.tolist().index(sel_val)
    #print idx, sel_val
    if idx != col:
        # change sign to keep the same det
        T[idx], T[col] = T[col].copy(), -T[idx]
        c[idx], c[col] = c[col], -c[idx]

    T[col] = dot(c,T)

    if c[col] < 0: # sign of det was changed, change it again
        T[1] *= -1

    if double_T:
        T /= 2.

    return T


def is_integer(a, epsilon=1e-7):
    "return true if numpy Float array consists off all integers"
    return (numpy.abs(a - numpy.round(a)) < epsilon).all()

def find_smallest_multiplier(a, max_n=1000):
    """return the smallest positive integer n such that matrix a multiplied
       by n is an integer matrix
    """
    for i in range(1, max_n):
        if is_integer(i*a):
            return i
    raise ValueError("Sorry, we can't make this matrix integer:\n%s" % a)

def find_smallest_real_multiplier(a, max_n=1000):
    """return the smallest positive real f such that matrix `a' multiplied
       by f is an integer matrix
    """
    # |the smallest non-zero element|
    m = min(abs(i) for i in a if abs(i) > 1e-9)
    for i in range(1, max_n):
        t = i / float(m)
        if is_integer(t * a):
            return t
    raise ValueError("Sorry, we can't make this matrix integer:\n%s" % a)

def scale_to_integers(v):
    return array(v * find_smallest_real_multiplier(v)).round().astype(int)

@transpose_3x3
def make_csl_from_0_lattice(T, n):
    if n < 0:
        T[0] *= -1
        n *= -1
    while True:
        m = [find_smallest_multiplier(T[i]) for i in (0,1,2)]
        prod = m[0] * m[1] * m[2]
        #print "prod", prod, n
        if prod <= n:
            for i in range(3):
                T[i] *= m[i]
            if prod < n:
                assert n % prod == 0
                T[0] *= n / prod
            break
        else:
            changed = False
            for i in range(3):
                for j in range(3):
                    if changed or i == j or m[i] == 1 or m[j] == 1:
                        continue
                    if m[i] <= m[j]:
                        a, b = i, j
                    else:
                        a, b = j, i
                    for k in plus_minus_gen(m[b]):
                        if find_smallest_multiplier(T[a] + k * T[b]) < m[a]:
                            T[a] += k * T[b]
                            changed = True
                            break
            assert changed, "Problem when making CSL from 0-lattice"
    assert is_integer(T)
    return T.round().astype(int)


def find_csl_matrix(sigma, R):
    """\
    Find matrix that determines the coincidence site lattice
    for cubic structures.
    Parameters:
        sigma: CSL sigma
        R: rotation matrix
        centering: "f" for f.c.c., "b" for b.c.c. and None for p.c.
    Return value:
        matrix, which column vectors are the unit vectors of the CSL.
    Based on H. Grimmer et al., Acta Cryst. (1974) A30, 197
    """

    S = _get_S()

    Rs = dot(dot(inv(S), inv(R)), S)
    #print "xxx",  inv(R)
    #print Rs
    found = False
    # searching for unimodular transformation that makes det(Tp) != 0
    for U in _get_unimodular_transformation():
        assert det(U) in (1, -1)
        Tp = identity(3) - dot(U, Rs)
        if abs(det(Tp)) > 1e-6:
            found = True
            print "Unimodular transformation used:\n%s" % U
            break
    if not found:
        print "Error. Try another unimodular transformation U to calculate T'"
        sys.exit(1)

    Xp = numpy.round(inv(Tp), 12)
    print "0-lattice:\n%s" % Xp
    n = round(sigma / det(Xp), 7)
    # n is an integral number of 0-lattice units per CSL unit
    print "det(X')=",det(Xp), "  n=", n
    csl = make_csl_from_0_lattice(Xp, n)
    assert is_integer(csl)
    csl = csl.round().astype(int)
    return beautify_matrix(csl)


def plus_minus_gen(n):
    for i in xrange(1, n):
        yield i
        yield -i

def zero_plus_minus_gen(n):
    yield 0
    for i in plus_minus_gen(n):
        yield i


@transpose_3x3
def find_orthorhombic_pbc(M):
    """\
     We don't change the last axis,
     because it was set properly in make_parallel_to_axis().

     Let M = [x, y, z], pbc = [x2, y2, z2],
     where M and pbc are matrices 3x3, and x, y, z, x2, y2, z2 are vectors.
     We search for multipliers b,c,d,e,f,g such that pbc is integer matrix,
         x2 = b x + d y + e z
         y2 = f x + c y + g z
         z2 = 0   + 0   + 1 z
     z2 = z/GCD(z), b != 0, c != 0,
     and max(|x2|, |y2|, |z2|) has the smallest value possible.
     (This description is a bit simplified, for details see the code)
     
     In matrix notation:       [b f 0]
                           M x [d c 0] = pbc  
                               [e g 1]        
    """
    # M is "half-integer" when using pc2fcc().
    # BTW I'm not sure if pc2fcc() is working properly.
    doubleM = not is_integer(M)
    if doubleM:
        M *= 2 # make it integer, will be /=2 at the end

    assert is_integer(M), M
    M = M.round().astype(int)

    # We are searching for solution by iteration over possible b,d,c values,
    # -max_multiplier < b,d,c < max_multiplier.
    # Increasing max_multiplier obviously slows down the program.
    max_multiplier = 27

    pbc = None
    max_sq = 0
    x, y, z = M

    x_ = x / gcd_array(x)
    y_ = y / gcd_array(y)
    z_ = z / gcd_array(z)
    #print "gcd_array", z, gcd_array(z)
    Mx = array([x, y_, z_])
    My = array([x_, y, z_])
    #print "mx",Mx
    #print "my",My

    z2 = z_ # previously: z2 = z
    mxz = dot(Mx, z2)
    myz = dot(My, z2)
    for b in plus_minus_gen(max_multiplier):
        for d in zero_plus_minus_gen(max_multiplier):
            e_ = - (mxz[0] * b + mxz[1] * d) / float(mxz[2])
            e = int(round(e_))
            if abs(e - e_) < 1e-7:
                x2 = dot([b,d,e], Mx)
                mxy = dot(My, x2)
                aa = array([[mxy[0],mxy[2]],
                            [myz[0],myz[2]]])
                bb = array([-mxy[1], -myz[1]])
                aa_invertible = (abs(det(aa)) > 1e-7)
                for c in plus_minus_gen(max_multiplier):
                    #  z2 . y2 == 0 and x2 . y2 == 0 =>
                    #       f * mxy[0] + g * mxy[2] == -c * mxy[1] 
                    #       f * myz[0] + g * myz[2] == -c * myz[1] 
                    if aa_invertible:
                        fg = solve(aa, c * bb)
                    else: # special case, i'm not sure if handled properly
                        for f in zero_plus_minus_gen(max_multiplier):
                            g_ = - (myz[0] * f + myz[1] * c) / float(myz[2])
                            g = int(round(g_))
                            if abs(g - g_) < 1e-7:
                                y2 = dot([f,c,g], My)
                                if inner(x2, y2) == 0:
                                    fg = array([f, g])
                                    break
                        else:
                            continue

                    if is_integer(fg) and (
                                  numpy.abs(fg) < max_multiplier - 0.5).all():
                        f, g = fg.round().astype(int)
                        y2 = dot([f,c,g], My)
                        max_sq_ = max(dot(x2,x2), dot(y2,y2), dot(z2,z2))
                        if pbc is None or max_sq_ < max_sq:
                            pbc = array([x2, y2, z2])
                            max_sq = max_sq_
    if pbc is None:
        print "No orthorhombic PBC found. (you may increase max_multiplier)"
        sys.exit()

    if doubleM:
        pbc /= 2.

    # we prefer determinant to be positive (negative one would later cause
    # inversion in addition to rotation)
    if det(pbc) < 0:
        pbc[0] = -pbc[0]

    # optionally swap x2 with y2
    id = identity(3)
    if (pbc[1] == id[0]).all() or (pbc[0] == -id[1]).all():
        pbc[0], pbc[1] = pbc[1], -pbc[0]
    elif (pbc[1] == -id[0]).all() or (pbc[0] == id[1]).all():
        pbc[0], pbc[1] = -pbc[1], pbc[0]

    return pbc


def find_type(type, Cp):
    for i,j,k in (0,0,1), (0,1,0), (1,0,0), (0,1,1), (1,0,1), (1,1,0), (1,1,1):
        if ((i * Cp[0] + j * Cp[1] + k * Cp[2]) % 2 == type).all():
            return [i,j,k]
    raise ValueError("find_type: %s not found" % type)

# see the paper by Grimmer, 1.3.1-1.3.3
@transpose_3x3
def pc2fcc(Cp):
    t1 = find_type([0,1,1], Cp)
    t2 = find_type([1,0,1], Cp)
    pos1 = t1.index(1)
    pos2 = t2.index(1)
    if pos2 == pos1:
        try:
            pos2 = t2.index(1, pos1+1)
        except ValueError:
            pos1 = t1.index(1, pos1+1)
    Z = identity(3)
    Z[pos1] = array(t1) / 2.
    Z[pos2] = array(t2) / 2.
    #print_matrix("Z (in pc2fcc)", Z.transpose())
    return dot(Z, Cp)

def print_list(hkl, max_angle, limit):
    print "[max. sigma: %s, max angle: %s deg.]" % (limit, max_angle)
    data = []
    for i in range(limit):
        tt = get_theta_m_n_list(hkl, i, verbose=False)
        for t in tt:
            theta, m, n = t
            if degrees(theta) <= max_angle:
                tup = (i, degrees(theta), m, n)
                data.append(tup)
                print "sigma=%3i    theta=%5.2f     m=%3i    n=%3i" % tup

    data.sort(key= lambda x: x[1])
    print " ============= Sorted by theta ================ "
    for i in data:
        print "sigma=%3i    theta=%5.2f     m=%3i    n=%3i" % i


def print_details(hkl, m, n):
    sigma = get_cubic_sigma(hkl, m, n)
    theta = get_cubic_theta(hkl, m, n)
    print "sigma=%d, theta=%.3f, m=%d, n=%d, axis=[%d,%d,%d]" % (
            sigma, degrees(theta), m, n, hkl[0], hkl[1], hkl[2])
    R = rodrigues(hkl, theta)
    print
    print "R * sigma =\n%s" % (R * sigma)
    C = find_csl_matrix(sigma, R)
    print "CSL primitive cell (det=%s):\n%s" % (det(C), C)
    ## optional, for FCC
    #C = pc2fcc(C)
    #C = beautify_matrix(C)
    #print_matrix("CSL cell for fcc:", C)

    Cp = make_parallel_to_axis(C, col=2, axis=hkl)
    if (Cp != C).any():
        print "after making z || %s:\n%s" % (hkl, Cp)
    pbc = find_orthorhombic_pbc(Cp)
    print_matrix("Minimal(?) orthorhombic PBC", pbc)


def main():
    # parse keyword options
    limit=1000
    max_angle=90
    for a in sys.argv[1:]:
        if '=' in a:
            key, value = a.split("=", 1)
            if key == "limit":
                limit = int(value)
            elif key == "max_angle":
                max_angle = float(value)
            else:
                raise KeyError("Unknown option: " + key)

    args = [a for a in sys.argv[1:] if '=' not in a]
    if len(args) < 1 or len(args) > 3:
        print usage_string
        return

    hkl = parse_miller(args[0])

    if len(args) == 1:
        print_list(hkl, max_angle=max_angle, limit=limit)
    elif len(args) == 2:
        sigma = int(args[1])
        thetas = get_theta_m_n_list(hkl, sigma, verbose=False)
        thetas.sort(key = lambda x: x[0])
        for theta, m, n in thetas:
            print "m=%2d  n=%2d %7.3f" % (m, n, degrees(theta))
        if not thetas:
            print "Not found."
    elif len(args) == 3:
        m = int(sys.argv[2])
        n = int(sys.argv[3])
        try:
            print_details(hkl, m, n)
        except KeyboardInterrupt:
            print " Interrupted. Exiting."



if __name__ == '__main__':
    main()


