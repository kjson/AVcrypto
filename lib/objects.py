
from sympy import *
from utils import * 
from random import randint 

"""

Classes for algebraic objects

"""

class AlgebraicSet(object):
	"""
	Algebraic set defined by polynomials over a field 
	Currently only finite fields of prime order are supported.
	"""

	def __init__(self, equations,field):
		if field.characteristic() <= 3:
			raise NotImplementedError("Only finite fields of prime order are currently supported.")

		symbols 		= set()
		polynomials 	= list()

		for f in equations:
			symbols = set(symbols | f.atoms(Symbol))
			polynomials.append(poly(f,symbols,domain=field))

		self.equations 	= equations
		self.field 		= field
		self.symbols 	= symbols 
		self.polynomials= polynomials

	def is_point(self,point):
		""" Returns true if point is zero set of defining functions """
		if len(point) != len(self.symbols):
			raise Exception("Given point is not it %s" % self.affine_space()) 

		for f in self.polynomials:
			if f.eval(point) != 0:
				return False 

		return True  

	def points(self):
		""" Returns points which satisfy all defining polynomials """
		return solve_poly_system(self.polynomials, self.field)

	def count_points(self):
		raise NotImplementedError("Still figuring out this http://ac.els-cdn.com/S0747717101904705/1-s2.0-S0747717101904705-main.pdf?_tid=2443a374-ccfa-11e4-beb9-00000aacb35e&acdnat=1426633725_4c513d8c3a848312863104ac79f4ce58")
		""" 
		Number of points satisfying all equations with the methods of 

		http://ac.els-cdn.com/S0747717101904705/1-s2.0-S0747717101904705-main.pdf?_tid=604b8ef6-ccee-11e4-9e66-00000aab0f01&acdnat=1426628671_1e73adbeb0d2deb89d89b55aa63774d1

		"""
		pass 


class Variety(AlgebraicSet):
	""" An irrecucible algebraic set"""
	def __init__(self, equations, field):
		AlgebraicSet.__init__(self,equations,field)

	def reducible(self):
		raise NotImplementedError("still figuring out proper criteria over finite fields")
		""" 
		Returns true if algebraic set can be written as the 
		product of two proper irreducible subsets 
		"""
		pass 


class AbelianVariety(Variety):
	""" A complete group variety with an operation on its points"""
	def __init__(self, equations, field, operation):
		Variety.__init__(self,equations,field)
		self.operation = operation

	def add(self,p1,p2):
		return self.operation(p1,p2)
		

class EllipticCurve(AbelianVariety):
	"""
	An abelian variety E of genus g = 1 and dimension d = 1 over a field of order p 

	Worklist

	- SEA
	- random point 
	- 
	"""
	def __init__(self,equations,field):

		def operation((x1,y1),(x2,y2)):

			p = field.characteristic()
			a = self.polynomials[0].coeffs()[1]

			if x1 != x2:  
				s = (y1 - y2) * modInv(x1 - x2,p) 

			if x1 == x2: 
				if y1 == -y2: 
					return (symbols('O'),symbols('O')) 
				else:
					s = (3*x1**2 + a) * modInv(2*y1,p) 

			x3 = s**2 - x1 - x2 
			y3 = - y1 + s * (x3 - x1)

			return (x3 % p,y3 % p)

		AbelianVariety.__init__(self,equations,field,operation)

		if len(self.polynomials[0].coeffs())!=4:
			raise NotImplementedError("Only elliptic curves with non-zero coeffiencts are supported.")

	def j_invariant(self):
		""" Computes the J-invarient of E """
		a = self.polynomials[0].coeffs()[1]
		b = self.polynomials[0].coeffs()[3]
		return 1728*(4*a**3/(4*a**3 + 27*b**2))

	def scalar_mult(self,num,point):
		""" Scalar multiplication of a point on E """
		if num == 0:
			return (symbols('O'),symbols('O'))
	 	if num == 1:
			return point
		if num % 2 == 0:
			return self.scalar_mult(num / 2, self.add(point,point))
		else:
			return self.add(point,self.scalar_mult(num - 1,point))

	def order(self):
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
		while m < 4*sqrt(q):
			if Elkies(l):
				ll = eigen_value()
				t = (ll + q/ll) % l 
				El.append((t,l))
			else:
				for t in possible_traces(l):
					A.append((t,l))

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

	def division_polynomial(self,n):
	    """ 
	    Calculate the nth division polynomial for E 
	    See http://en.wikipedia.org/wiki/Division_polynomials#Definition
	    for definition of the nth division polynomial.
	    """

	    # Indeterminants 
	    x,y = symbols('x,y')

	    # Alias function for simplicity 
	    f = division_polynomial

	    # Get coefficients
	    a = self.polynomials[0].coeffs()[1]
	    b = self.polynomials[0].coeffs()[3]

		# index n must be possitive
	    assert n >= 0 

	    # Base cases for recursion 
	    if n==0 or n==1:
	        return n 
	    if n==2:
	        return 2*y 

	    # Haven't worked out if we need these yet... probably wont 
	    if n==3:
	        return 3*x**4 + 6*a*x**2 + 12*b*x - a**2 
	    if n==4:
	        return 4*y*(x**6 + 5*a*x**4 + 20*b*x**3 - 5*a**2*x**2 - 4*a*b*x - 8*b**2 -a**3)

	    # m > 1
	    if n%2==1:
	        m = n/2 
	        return f(a,b,m+2)*f(a,b,m)**3 - f(a,b,m-1)*f(a,b,m+1)**3 
	    # m > 2 
	    if n%2 ==0:
	        m = n/2 
	        return (f(a,b,m) / 2*y)*(f(a,b,m+2)*f(a,b,m-1)**2 - f(a,b,m-2)*f(a,b,m+1)**2)

		

class HyperEllipticCurve(AbelianVariety):
	""" Class for HyperEllipticCurves """
	def __init__(self, equations,field):

		def operation((a1,b1),(a2,b2)):	

			f = poly(equations[0],domain=field)

			# TEmporraily set 
			g = 3

			# Composition 

			d1, e1, c1 	= gcdex(a1,a2)
			d, e2, c2 	= gcdex(d1, b1 + b2)

			a3 = (a1*a2) / d**2 

			# b_3 = reduce(((c2*c1*a1 + c2*e1*a2 + e2*(b1*b2 + f)) / d),[a_3])[1]
			print reduce(((c2*c1*a1 + c2*e1*a2 + e2*(b1*b2 + f)) / d),[a3])

			# Reduction 

			# u2 = (f - v1**2) / u1 
			# v2 = reduced((-v1), [u2])[1]

			# while degree(u2) > g:
			# 	u2 = (f - v2**2) / u2
			# 	v2 = reduced((-v2), [u2])[1]

			# if coeffs(u2)[0] != 1:
			# 	u2 = u2/ coeffs(u2)[0]

			# u2 = simplify(u2)
			# v2 = simplify(v1)

			# return (a3,b3)

		AbelianVariety.__init__(self,equations,field,operation)







