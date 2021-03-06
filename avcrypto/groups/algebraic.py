import sympy as sp
import random as rn
import math
from Util.number_theory import extended_gcd, mod_inv, is_probable_prime, \
     lenstra, legendre, chinese_remainder_theorem


class Group(object):
    """ Abstract abelian finite group """
    def __init__(self, elements,operation):
        self.elements = points
        self.operation = operation
        
    def __mul__(self,other):
        """ Cartesian product of two groups. Operation is the 
        canonically cartesian operation """
        def operation((x1,y1),(x2,y2)):
            return (self.operation(x1,x2),other.operation(y1,y2))
        return Group((self.points,other.points),operation)

    def is_point(self,point):
        """ Checks if point is in group """   
        if point in self.points:
            return True
        else:
            return False

    def scalar_multiplication(self,scalar,point):
        """ naive scalar multiplication of points in group.
        Hopefully a child class of Group has a better implentation. """
        for i in xrange(0,scalar):
            point = point + point
        return point

class Point(object):
    """ Point in an abelian finite group """
    def __init__(self,point,group):
        assert group.is_point(point)
        self.group = group
        self.point = point
        
    def __add__(self,other):
        """ overload + operator with group operation """
        return self.group.operation(self,other)
  
    def __mul__(scalar,self):
        """ Scalar multiplication of a point """  
        return self.group.scalar_multiplication(scalar,self)    


class AbelianVariety(Group):
    """
    A complete algebraic abelian group defined as the zero locus of
    an ideal of polynomials over a finite field of prime order.
    Let A represent an Abelian Variety in all the doc-strings.
    """
    def __init__(self,equations,field,operation):

        if field.characteristic() <= 3: # Characteristic != 2,3
            raise NotImplementedError("Only finite fields of prime order are \
                supported.")

        symbols = set()
        polynomials = list() # Force sympy Poly class for equations
        for f in equations:
            symbols = set(symbols | f.atoms(sp.Symbol))
            polynomials.append(sp.poly(f, symbols, domain=field))

        self.equations = equations
        self.field = field
        self.operation = operation
        self.polynomials = polynomials
        self.symbols = symbols


        # Group.__init__(self,sp.solve_poly_system(self.polynomials,
        # self.field),operation)
        Group.__init__(self,[],operation)


    def points(self):
        """ Finds all points which satisfy eqautions of A """
        return sp.solve_poly_system(self.polynomials, self.field)

    def is_point(self, point):
        """ Verifies that given point belongs to A """
        if len(point) != len(self.symbols):
            raise Exception("Point must be of proper tuple size")
        for f in self.polynomials:
            if f.eval(point) != 0:
                return False
        return True

    def order(self):
        """ Returns number of points in A """

    def irreducible(self):
        """ Determines if A is irreducible """

    def add(self, p1, p2):
        """ Adds two points on A using defined operation """
        return self.operation(p1, p2)


class EllipticCurve(AbelianVariety):
    """ An elliptic curve E defined by an equation y**2 = x**3 + a*x + b """
    def __init__(self, equations, field):

        def operation((x1, y1), (x2, y2)):
            """ Addition of two points on E """
            p = field.characteristic()
            a = self.polynomials[0].coeffs()[1]
            if x1 != x2:
                s = (y1 - y2) * mod_inv(x1 - x2, p)
            if x1 == x2:
                if y1 == -y2:
                    return (sp.symbols('O'), sp.symbols('O'))
                else:
                    s = (3 * x1 ** 2 + a) * mod_inv(2 * y1, p)
            x3 = s ** 2 - x1 - x2
            y3 = - y1 + s * (x1 - x3)
            return (x3 % p, y3 % p)

        if type(equations) != 'list':
            equations = [equations]

        AbelianVariety.__init__(self, equations, field, operation)

        a = self.a = self.polynomials[0].coeffs()[1]
        b = self.b = self.polynomials[0].coeffs()[3]
        p = self.p = self.field.characteristic()

        # Check both a,b are not zero
        if len(self.polynomials[0].coeffs()) != 4:
            raise NotImplementedError("Only elliptic curves with non-zero \
                coeffiencts are supported.")

        # Curve cannot be singular
        if (4*a**3 + 27*b**2) % p == 0:
            raise Exception("Curve cannot be singular")

        # Only one equaions may be given
        if len(self.polynomials) > 1:
            raise Exception("Elliptic curve must be given by equations \
                y**2 = x**3 + a*x + b ")


    def scalar_mult(self, scalar, point):
        """ Scalar multiplication of a point on E """
        if scalar == 0:
            return (sp.symbols('O'), sp.symbols('O'))
        elif scalar == 1:
            return point
        elif scalar % 2 == 0:
            return self.scalar_mult(scalar / 2, self.add(point, point))
        else:
            return self.add(point, self.scalar_mult(scalar - 1, point))


    def random_point(self):
        """ Finds random a point (x,y) on E """
        a,b,p = self.a,self.b,self.p
        x = rn.randrange(1, p)

        while legendre(x ** 3 + a * x + b, p) != 1:
            x = rn.randrange(1, p)

        Z = x ** 3 + a * x + b

        if p % 4 == 3:
            return (x,(Z ** ((p + 1) / 4)) % p)
        else:
            r = rn.randrange(1, p)

            while legendre(r ** 2 - 4 * Z, p) != -1:
                r = rn.randrange(1, p)

            d = (r ** 2 - 4 * Z) % p
            alpha = (r + math.sqrt(d)) / 2
            beta  = (alpha ** ((p + 1) / 2)).expand()
            beta = str(beta)
            for string in re.split("([+-/])",beta.replace(" ","")):
                if "sqrt" in string:
                    for number in re.split("([*])",string):
                        if not "sqrt" in number:
                            if not "*" in number:
                                return ((int(number)/2) %p,x)

    def order(self):
        """ Order of E over F_p"""
        p = self.p
        if p.bit_length() < 20:
            return self.schoof()
        elif 20 <= p.bit_length() <= 100:
            return self.schoof()
        else:
            return self.schoof_elkies_atkin()

    def lenstra(self):
        """ Lenstra's simple point counting algorithm for elliptic curves"""
        a,b,p = self.a,self.b,self.p
        number_of_points = 1 + p
        for x in xrange(0,p):
            number_of_points+= legendre(x**3 + a*x + b, p)
        return number_of_points

    def schoof(self):
        """
        schoof's algorithm for counting the number of points on an elliptic
        curve over F_p

        The goal is to compute the trace of the frobenious endomorphism t. 
        Then p + 1 + t is the order of the group. 

        Computing t is done by computing t mod l_1,l_2,...,l_s for a 
        sufficient number of primes l_1,l_2,...,l_s such that their 
        product is greater than 4*sqrt(p).

        Then use the chinese remainder theorem to recover t. let q_li and t_li
        be q and t mod l_i for i= 1,...,s. Then the equation 
            
                (x**2p,y**2p) + q_li(x,y) = t_l(x**p,y**p) mod l_i
        
        is used to calculate the t_li's.  

        """

        a,b,p = self.a,self.b,self.p
        x, y = sp.symbols('x,y')
        list_of_primes=list_of_congruences=[]

        # Build list of odd primes satisfying condition's above
        i=3
        product = 2
        while product < 4*math.sqrt(p):
            if is_probable_prime(i):
                list_of_primes.append(i)
                product*=i
            i+=2

        # Special case to determine t mod 2. 
        if roots_in_F_q(p, f):
            list_of_congruences.append((2, 1))
        else:
            list_of_congruences.append((2, 0))

        # alias E.dpoly()
        f = E.dpoly

        # Build list of congruences
        for l in list_of_primes:
            p_l = p % l 
            x_pl,y_pl = x - f(pl-1)*f(pl+1) / f(pl)**2, f(2*pl) / f(pl)**4

            # We skip the case that (x**2p,y**2p) == +- q_pl(x,y)
            if (x**(2*p), y**(2*p)) == (x_pl, y_pl):
                continue
            elif (x**(2*p), y**(2*p)) == (-x_pl, -y_pl):
                continue

            # Calculate (x',y') = (x**2p,y**2p) + q_pl(x,y) mod phi_l
            x_prime, y_prime = self.symbolic_add((x**(2*p), y**(2*p)), (x_pl, y_pl))
            x_prime, y_prime = sp.reduced(x_prime,f(pl)),sp.reduced(y_prime,f(pl))

            # Look for t_l such that (x',y') = t_l(x**p,y**p) mod phi(l) 
            for t_l in xrange(1,(l-1)/2):
                found_tl = False
                x_tl,y_tl = self.symbolic_scalar(t_l,(x**p,y**p))
                if sp.reduced(x_prime - x_tl,f(pl)) == 0:
                    found_tl = True
                    if sp.reduced((y_prime - y_tl)/y, f(pl)) == 0:
                        list_of_congruences.append((l,t_l))
                    else: 
                        list_of_congruences.append((l,-t_l))

        return chinese_remainder_theorem(list_of_congruences)


    def dpoly(self, n):
        """ Calculate the nth division polynomial for E """
        x, y = sp.symbols('x,y')
        a,b = self.a,self.b
        assert n >= 0

        if n == 0 or n == 1:
            return sp.poly(n, gens=x)
        elif n == 2:
            return sp.poly(2 * y)
        elif n == 3:
            return sp.poly(3 * x ** 4 + 6 * a * x ** 2 + 12 * b * x - a ** 2)
        elif n == 4:
            return sp.poly(4 * y * (x ** 6 + 5 * a * x ** 4 \
                                   + 20 * b * x ** 3 - 5 * a ** 2 * x ** 2 - \
                                   4 * a * b * x - 8 * b ** 2 - a ** 3))
        elif n % 2 == 1:
            m = n / 2
            return self.dpoly(m + 2) * self.dpoly(m) ** 3 \
                   - self.dpoly(m - 1) * self.dpoly(m + 1) ** 3
        elif n % 2 == 0:
            m = n / 2
            return (self.dpoly(m) / 2 * y) \
                   * (self.dpoly(m + 2) * self.dpoly(m - 1) ** 2 \
                      - self.dpoly(m - 2) * self.dpoly(m + 1) ** 2)


    def symbolic_scalar(self, scalar, (x, y)):
        """ Symbolic scalar multiplication  of (x,y) """
        if scalar == 0:
            return (sp.symbols('O'), sp.symbols('O'))
        elif scalar == 1:
            return (x, y)
        elif scalar % 2 == 0:
            return self.symbolic_scalar(scalar / 2, self.symbolic_add((x, y), (x, y)))
        else:
            return self.symbolic_add((x, y), self.symbolic_scalar(scalar - 1, (x, y)))


    def symbolic_add(self, (x1, y1), (x2, y2)):
        """ Adds (x1,y1) + (x2,y2) symbolically """
        a = self.a
        if x1 != x2:
            s = (y1 - y2) / (x1 - x2)
        else:
            if y1 == -y2:
                return (sp.symbols('O'), sp.symbols('O'))
            else:
                s = (3 * x1 ** 2 + a) / (2 * y1)
        x3 = s ** 2 - x1 - x2
        y3 = - y1 + s * (x1 - x3)
        return (x3, y3)


    def sea(self):
        """" Schoof-Elkies-Atkin Algorithm for computing the number of points
        in E. This file uses modular polynomials and isogenies volcanos.

        Info found here ./.res/Implementing the Schoof-Elkies-Atkin Algorithm with NTL.pdf """

        a,b,p = self.a,self.b,self.p
        list_of_primes = Alkies = Elkies = []
        x,y = sp.symbols('x,y')


        # 1 Figure out first congruence l = 2
        if sp.gcd(x**p -x,x**3+a*x +b)!=1:
            Elkies.append((0,2))
        else:
            Elkies.append((1,2))

        # 3 Create list of primes
        i=3
        product = 2
        while product < 4*math.sqrt(p):
            if is_probable_prime(i):
                list_of_primes.append(i)
                product*=i
            i+=2


        # 3
        for l in list_of_primes:
            # a) Computer lth modular polynomial phi

            # a) evaluate phi at y = j (j is j-invariant of E), call this phi(j)
            # b) compute X**p mod phi(j)
            # c) Compute gcd(phi(j),x**P - x)
            # c) decide if l is an elkies or atkins prime
            if elkies == True:
                # d) i compute j(E) as the root of phi_lj in F_p
                # d)
                # d)
                # d)
                # d)
                pass 
            else:
                pass 





        t = match_sort(Elkies + Atkins)
        return p + 1 + t


    def order_frobenious(self,phi_lj):
        pass







class HyperEllipticCurve(Group):
    """
    A hyperelliptic curve defined by an equation 0 = y**2 - f(x)
    where deg(f) > 4, with no repeating roots over a finite field.
    Throughout these docstrings, the capital letter H will be used
    to denote a HyperEllipticCurve
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
            e1, e2, d1 = sp.gcdex(u1, u2)
            c1, c2, d = sp.gcdex(d1, v1 + v2)
            s1, s2, s3 = c1 * e1, c1 * e2, c2
            u = u1 * u2 / d ** 2
            v = sp.rem((s1 * u1 * v2 + s2 * u2 * v1 + s3 * (v1 * v2 + f)) / d, u)

            while sp.degree(u) > g:
                u = (f - v ** 2) / u
                v = sp.rem(-v, u)

            if u.coeffs()[len(u.coeffs()) - 1] != 1:
                u = u / u.coeffs()[len(u.coeffs()) - 1]

            return (u, v)

        def elements():
            pass

        Group.__init__(self,elements,operation)

    def order(self):
        """ Returns the number of points on H """
