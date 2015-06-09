import sympy as S
from groups import Group
from util import *

class HyperEllipticCurve(Group):
    """
    A hyperelliptic curve defined by an equation 0 = y**2 - f(x)
    where deg(f) > 4, with no repeating roots over a finite field.
    """

    def __init__(self, equations, field):

        def genus(self):
            """ The arithmetic genus of H """
            f = equations[0].subs(y ** 2, 0)
            if degree(f) == 3:
                return 1
            else:
                return 2

        self.genus = self.genus()

        def operation((u1, v1), (u2, v2)):
            """ Addition of two divisors using cantor's algorithm """

            f = equations[0].subs(y ** 2, 0)
            g = self.genus
            e1, e2, d1 = S.gcdex(u1, u2)
            c1, c2, d = S.gcdex(d1, v1 + v2)
            s1, s2, s3 = c1 * e1, c1 * e2, c2
            u = u1 * u2 / d ** 2
            v = S.rem((s1 * u1 * v2 + s2 * u2 * v1 + s3 * (v1 * v2 + f)) / d, u)

            while S.degree(u) > g:
                u = (f - v ** 2) / u
                v = S.rem(-v, u)

            if u.coeffs()[len(u.coeffs()) - 1] != 1:
                u = u / u.coeffs()[len(u.coeffs()) - 1]

            return (u, v)

        def elements():
            pass

        Group.__init__(self,elements,operation)

    def order(self):
        """ Returns the number of points in H using BLANK algorithm"""
