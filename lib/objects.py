from sympy import *
from utils import modInv
from random import randint 

class AlgebraicSet(object):
	"""

	Algebraic set defined as the zero locus of given

	--polynomials
	--field 

	Currently only finite fields or prime order are supported.
	
	"""

	def __init__(self, equations,field):
		super(AlgebraicSet, self).__init__()
		self.equations = equations
		self.field = field
		
		if field.characteristic() <= 3:
			raise NotImplementedError("Only finite fields of prime order are supported.")

		symbols = set()
		polynomials = list()

		for f in equations:
			symbols = set(symbols | f.atoms(Symbol))
			polynomials.append(poly(f,symbols,domain=field))

		self.symbols = symbols 
		self.polynomials = polynomials

	def is_point(self,point):
		""" Returns true if point is zero set of defining functions """
		if len(point) != len(self.symbols):
			raise Exception("Given point is not it %s" % self.affine_space()) 

		for f in self.polynomials:
			if f.eval(point) != 0:
				return False 
		return True  

	def set_polynomials(self,equations):

		polynomials = list()
		symbols 	= set()

		for f in equations:
			symbols = set(symbols | f.atoms(Symbol))
			polynomials.append(poly(f,symbols,domain=self.field))

		self.symbols = symbols
		self.polynomials = polynomials


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
			y3 = y1 + s * (x3 - x1)

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

class HyperEllipticCurve(AbelianVariety):
	"""docstring for HyperEllipticCurve"""
	def __init__(self, equations,field):

		def operation((a1,b1),(a2,b2)):	
			d1 	= gcd(a1,a2)
			d 	= gcd(d1, b1 + b2)

			# handbook pg. 308

			

		AbelianVariety.__init__(self,equations,field,operation)


	def number_of_points(self):
		pass









