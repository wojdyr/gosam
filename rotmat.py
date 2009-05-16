# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2
"""\
Mathematical and matrix-related utilities.
Most notably Rodrigues rotation formula.
"""

from math import cos, sin, acos, asin, sqrt, pi, radians, degrees
from numpy import array, dot, inner, identity, linalg, sign

def rodrigues(a, angle, verbose=False):
    "use Rodrigues' rotation formula to get rotation matrix"
    a = array(a, dtype=float)
    a /= sqrt(inner(a, a)) # make unit vector
    #assert abs(sin_angle - sin(acos(cos_angle))) < 1e-6
    if verbose:
        print "rotation angle:", degrees(angle)
        print "rotation axis:", a
    omega = array([[   0., -a[2],  a[1]],
                   [ a[2],    0., -a[0]],
                   [-a[1],  a[0],    0.]])
    rm = (identity(3) + omega * sin(angle)
                            + dot(omega, omega) * (1 - cos(angle)))
    if verbose:
        print "rotation matrix:", rm
    return rm


def print_matrix(text, M):
    print "%s: (det=%s):\n%s" % (text, linalg.det(M), M)


def round_to_multiplicity(m, val):
    "round val to the nearest multiplicity of m, but avoid zero"
    return (round(float(val) / m) or 1) * m

def is_diagonal(m):
    m = array(m)
    assert m.shape == (3,3)
    return (sign(m) == identity(3)).all()


class StdDev:
    def __init__(self):
        self.n = 0
        self.mean = 0.
        self.S = 0.

    def __str__(self):
        return "%s +- %s" % (self.mean, self.get_stddev())

    def get_variance(self):
        return self.S / (self.n - 1)

    def get_stddev(self):
        return sqrt(self.get_variance())

    def add_x(self, x):
        self.n += 1
        delta = x - self.mean
        self.mean += delta / self.n
        self.S += delta * (x - self.mean)


