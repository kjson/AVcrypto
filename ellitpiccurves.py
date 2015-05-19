from abelianvarieties import AbelianVariety
from util import *
import sympy as S  
import re 

class EllipticCurve(AbelianVariety):
    """ A elliptic curve E defined by an equation y**2 = x**3 + a*x + b """

    def __init__(self, equations, field):

        def operation((x1, y1), (x2, y2)):
            """ Addition of two points on E """
            p = field.characteristic()
            a = self.polynomials[0].coeffs()[1]

            if x1 != x2:
                s = (y1 - y2) * modInv(x1 - x2, p)

            if x1 == x2:
                if y1 == -y2:
                    return (S.symbols('O'), S.symbols('O'))
                else:
                    s = (3 * x1 ** 2 + a) * modInv(2 * y1, p)

            x3 = s ** 2 - x1 - x2
            y3 = - y1 + s * (x3 - x1)

            return (x3 % p, y3 % p)

        AbelianVariety.__init__(self, equations, field, operation)

        # Check both a,b are not zero
        if len(self.polynomials[0].coeffs()) != 4:
            raise NotImplementedError("Only elliptic curves with non-zero coeffiencts are supported.")

    def j_invariant(self):
        """ Computes the J-invarient of E """
        a = self.polynomials[0].coeffs()[1]
        b = self.polynomials[0].coeffs()[3]
        return 1728 * (4 * a ** 3 / (4 * a ** 3 + 27 * b ** 2))

    def random_point(self):
        """ Finds random point on E"""

        a = self.polynomials[0].coeffs()[1]
        b = self.polynomials[0].coeffs()[3]
        p = self.field.characteristic()

        x = random.randrange(1, p)

        while legendre(x ** 3 + a * x + b, p) != 1:
            x = random.randrange(1, p)

        if p % 4 == 3:
            return (((x ** 3 + a * x + b) ** ((p + 1) / 4)) % p, x)
        else:
            r = random.randrange(1, p)
            while legendre(r ** 2 - 4 * (x ** 3 + a * x + b), p) != -1:
                r = random.randrange(1, p)

            d = r ** 2 - 4 * (x ** 3 + a * x + b) % p 
            alpha = (r + sqrt(d)) / 2
            beta  = (alpha ** ((p + 1) / 2)).expand()

            # clean this up: You just want to int in front of sqrt(d)
            beta = str(beta)
            for string in re.split("([+-/])",beta.replace(" ","")):
                if "sqrt" in string:
                    for number in re.split("([*])",string):
                        if not "sqrt" in number:
                            if not "*" in number: 
                                return ((int(number)/2) %p,x)

    def scalar_mult(self, num, point):
        """ Scalar multiplication of a point on E """
        if num == 0:
            return (S.symbols('O'), S.symbols('O'))
        if num == 1:
            return point
        if num % 2 == 0:
            return self.scalar_mult(num / 2, self.add(point, point))
        else:
            return self.add(point, self.scalar_mult(num - 1, point))


    def order(self):
        """ Order of the group of points """
        p = self.field.characteristic()
        if bitlength(p) < 20:
            return self.lenstra()
        elif 20 <= bitlength(p) <= 100:
            return self.schoof()
        else: 
            return self.sea()


    def lenstra(self):
        """ Lenstra's simple point counting method """
        p = self.field.characteristic()
        a = self.polynomials[0].coeffs()[1]
        b = self.polynomials[0].coeffs()[3]
        number_of_points = 1 + p 
        for x in xrange(1,p):
            number_of_points += legendre((x**3 + a*x + b) % p,p)
        return number_of_points


    def schoof(self):
        """ Rene Schoofs point counting algorithm """

        # Values associated to E
        p = self.field.characteristic()
        a = self.polynomials[0].coeffs()[1]
        b = self.polynomials[0].coeffs()[3]
        x, y = S.symbols('x,y')

        # Build set of primes S
        listOfPrimes = [2, 3]
        i = 3
        s = 6
        limit = 4 * math.sqrt(p)
        while s < limit:
            if is_probable_prime(i + 2):
                listOfPrimes.append(i + 2)
                s *= i + 2
            i += 1

        # Build list of congruences
        listOfCongruences = list()

        if S.gcd(x ** p - x, x ** 3 + a * x + b) != 1:
            listOfCongruences.append((2, 0))
        else:
            listOfCongruences.append((2, 1))

        for l in listOfPrimes[1:]:
            div_poly_l = E.dpoly(l).subs(y ** 2, x ** 3 - a * x + b)
            pl = p % l
            x_pl = (x - E.dpoly(pl - 1) * E.dpoly(pl + 1) / (E.dpoly(pl)) ** 2).subs(y ** 2, x ** 3 - a * x + b)
            y_pl = (E.dpoly(2 * pl) / E.dpoly(pl) ** 4).subs(y ** 2, x ** 3 - a * x + b)
            theta_y = y_pl / y
            # Problem polynomial has massive degree
            slope = ((x ** 3 + a + x + b) * ((x ** 3 + a * x + b) ** (((p ** 2) - 1) / 2) - theta_y ) ** 2) / (
                x ** (2 * p) - x_pl)
            x_prime = slope ** 2 - x ** (2 * p) - x_pl
            y_prime = - (x ** 3 + a + x + b) ** p + slope * (x_prime - x ** (2 * p))

            for t in xrange(1, (l - 1) / 2 + 2):
                rhs = E.symbolic_scalar(t, (x, y))[0]
                if S.div(x_prime - rhs, div_poly_l)[1] == 0:
                    if S.div((y_prime - E.dpoly(2 * t) / 2 * E.dpoly(t) ** 4) / y, div_poly_l)[1] == 0:
                        listOfCongruences.append((l, t))
                    else:
                        listOfCongruences.append((l, -t))

            
        print listOfCongruences

    def symbolic_scalar(self, num, (x, y)):
        """ Symbolic scalar multiplication of points """
        if num == 0:
            return (S.symbols('O'), S.symbols('O'))
        if num == 1:
            return (x, y)
        if num % 2 == 0:
            return self.symbolic_scalar(num / 2, self.symbolic_add((x, y), (x, y)))
        else:
            return self.symbolic_add((x, y), self.symbolic_scalar(num - 1, (x, y)))

    def symbolic_add(self, (x1, y1), (x2, y2)):
        """ Adds (x,y) to itself symbolically """
        a = self.polynomials[0].coeffs()[1]
        if x1 != x2:
            s = (y1 - y2) / (x1 - x2)
        else:
            if y1 == -y2:
                return (S.symbols('O'), S.symbols('O'))
            else:
                s = (3 * x1 ** 2 + a) / (2 * y1)
        x3 = s ** 2 - x1 - x2
        y3 = - y1 + s * (x3 - x1)

        return (x3, y3)

    def sea(self):
        raise NotImplementedError("See worklist.")
        """
        Determine the number of points of E over F
        using the Schoof-Elkies-Atkin Algorithm.

        See reference A Brief Introduction to the Scoof-Elkies-Atkins(SEA) Algorithm.pdf

        WorkList

        -Supersingular checks  and p > 2l + 2
        -Elkies prime
        -Eigenvalues of map of frobenious
        -calc lth modular polynomial
        -fast factorization of univariate polynomials over finite fields see gnfs
        """

        # All variables used in algorithm
        m, l, A, El, q = 1, 2, list(), list(), self.field.characteristic

        # Only works for smooth curves
        assert self.j_invariant() != 0 or self.j_invariant() != 1728

        # Find possible traces
        while m < 4 * sqrt(q):
            if Elkies(l):
                ll = eigen_value()
                t = (ll + q / ll) % l
                El.append((t, l))
            else:
                for t in possible_traces(l):
                    A.append((t, l))

            m *= l
            l = next_prime(l)

        # return number of points on E.
        return q + 1 - t

        # These are really specific to order finding.
        def Elkies(l):
            """
            Determine if l is a Elkies prime or not

            """
            pass

        def eigen_value():
            """ Find Igenvalue
            See report/references/Fast algorithms for computing the eigenvalue in the Schoof-Elkies-Atkin algorithm.ps
            """
            pass

        def possible_traces(l):
            """ Find all traces of mod l """
            pass

    def dpoly(self, n):
        """
        Calculate the nth division polynomial for E
        See http://en.wikipedia.org/wiki/dpolys#Definition
        for definition of the nth division polynomial.
        """

        # Indeterminants
        x, y = S.symbols('x,y')

        # Get coefficients
        a = self.polynomials[0].coeffs()[1]
        b = self.polynomials[0].coeffs()[3]

        # index n must be possitive
        assert n >= 0

        # Base cases for recursion
        if n == 0 or n == 1:
            return S.poly(n, gens=x)
        if n == 2:
            return S.poly(2 * y)
            # Haven't worked out if we need these yet... probably wont
        if n == 3:
            return S.poly(3 * x ** 4 + 6 * a * x ** 2 + 12 * b * x - a ** 2)
        if n == 4:
            return S.poly(4 * y * (x ** 6 + 5 * a * x ** 4 \
                                   + 20 * b * x ** 3 - 5 * a ** 2 * x ** 2 - \
                                   4 * a * b * x - 8 * b ** 2 - a ** 3))
        if n % 2 == 1:
            m = n / 2
            return self.dpoly(m + 2) * self.dpoly(m) ** 3 \
                   - self.dpoly(m - 1) * self.dpoly(m + 1) ** 3
        if n % 2 == 0:
            m = n / 2
            return (self.dpoly(m) / 2 * y) \
                   * (self.dpoly(m + 2) * self.dpoly(m - 1) ** 2 \
                      - self.dpoly(m - 2) * self.dpoly(m + 1) ** 2)



x, y = S.symbols('x,y')
F = S.FiniteField(41)
f =  y ** 2 - x ** 3 - x - 1
E = EllipticCurve([f], F)
print E.lenstra()





