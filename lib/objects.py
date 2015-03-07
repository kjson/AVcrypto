from sympy import *
from utils import * 
from random import randint 

class AlgebraicSet(object):

	"""
	Algebraic set defined by polynomials over a field 
	Currently only finite fields of prime order are supported.
	"""

	def __init__(self, equations,field):
		if field.characteristic() <= 3:
			raise NotImplementedError("Only finite fields of prime order are supported.")

		symbols 	= set()
		polynomials = list()

		for f in equations:
			symbols = set(symbols | f.atoms(Symbol))
			polynomials.append(poly(f,symbols,domain=field))

		self.equations 		= equations
		self.field 			= field
		self.symbols 		= symbols 
		self.polynomials 	= polynomials

	def is_point(self,point):
		""" Returns true if point is zero set of defining functions """
		if len(point) != len(self.symbols):
			raise Exception("Given point is not it %s" % self.affine_space()) 

		for f in self.polynomials:
			if f.eval(point) != 0:
				return False 
		return True  

	def points(self):
		return solve_poly_system(self.polynomials)

	# def set_polynomials(self,equations):

	# 	polynomials = list()
	# 	symbols 	= set()

	# 	for f in equations:
	# 		symbols = set(symbols | f.atoms(Symbol))
	# 		polynomials.append(poly(f,symbols,domain=self.field))

	# 	self.symbols = symbols
	# 	self.polynomials = polynomials


class Variety(AlgebraicSet):

	""" An irrecucible algebraic set"""

	def __init__(self, equations, field):
		AlgebraicSet.__init__(self,equations,field)



class AbelianVariety(Variety):

	""" A complete group variety with an operation on its points"""

	def __init__(self, equations, field, operation):
		Variety.__init__(self,equations,field)
		self.operation = operation

	def add(self,p1,p2):
		return self.operation(p1,p2)
		

class EllipticCurve(AbelianVariety):
	"""An abelian variety of genus 1 and dimension 1"""
	def __init__(self,equations,field):

		def operation((x1,y1),(x2,y2)):

			p = field.characteristic()
			a = self.polynomials[0].coeffs()[2]

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

	def get_points(self):
		p = self.field.characteristic()
		list_of_points = list()

		for x in xrange(0,p):
			for y in xrange(0,p):
				if self.is_point((x,y)):
					list_of_points.append((x,y))

		return list_of_points

	def random_point(self):
		p = self.field.characteristic()
		x = randint(0,p)
		y = randint(0,p)

		while not self.is_point((x,y)):
			x = randint(0,p)
			y = randint(0,p)

		return (x,y)

	def scalar_mult(self,num,point):
		if num == 0:
			return (symbols('O'),symbols('O'))
	 	if num == 1:
			return point
		if num % 2 == 0:
			return self.scalar_mult(num / 2, self.add(point,point))
		else:
			return self.add(point,self.scalar_mult(num - 1,point))


	def order(self):
		""" Schoof's algorithm to count number of points on elliptic curve """

		p = self.field.characteristic()

		N = 4*sqrt(p) + 1 

		L = easy_factor(N)
		traces = list()

		for i in xrange(0,len(L)-1):
			traces[i] = trace_frob(p,L[i])

		return p + 1 - chinese_remainder_theorem(L,traces,len(L))



		# use  Schoof-Elkies-Atkin Algorithm 



class HyperEllipticCurve(AbelianVariety):
	""" HyperEllipticCurves """
	def __init__(self, equations,field):

		def operation((a1,b1),(a2,b2)):	

			f = poly(equations[0],domain=field)

			# TEmporraily set 
			g = 3

			# Composition 

			d1, e1, e2 	= gcdex(a1,a2)
			d2, c1, c2 	= gcdex(d1, b1 + b2)

			s1 = c1*e1
			s2 = c1*e2
			s3 = c2 

			u1 = (a1*a2) / d2**2 
			v1 = reduced(s1*a1*b2 + s1*a2*b1 + s3*b1*b2 + s3*f, [u1])[1]

			u1 = simplify(u1)
			v1 = simplify(v1)

			# Reduction 

			u2 = (f - v1**2) / u1 
			v2 = reduced((-v1), [u2])[1]

			while degree(u2) > g:
				u2 = (f - v2**2) / u2
				v2 = reduced((-v2), [u2])[1]

			if coeffs(u2)[0] != 1:
				u2 = u2/ coeffs(u2)[0]

			u2 = simplify(u2)
			v2 = simplify(v1)

			return (u2,v2)

		AbelianVariety.__init__(self,equations,field,operation)


	def number_of_points(self):
		pass

	

x,y = symbols('x,y')


f = y**2 - x**3 - x - 1 
g = x*y + y**2 + x 

F = FiniteField(43)

E = EllipticCurve([f],F)

print E.order()




