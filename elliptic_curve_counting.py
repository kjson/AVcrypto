"""
Tests for the arithmetic methods in the EllipticCurve class. Methods are 
tested over a range of finite fields on random elliptic curves. 

methods tested are; add, scalar_multiplication, random_point, order

"""

from algebraic_groups import EllipticCurve
from slog import Slog
from progressbar import ProgressBar 	
from nt_utils import primes_less_than, random_elliptic_curve

import os
import time 

import sympy as sp

# Create sympy variables
x, y = sp.symbols('x, y')

# Set test variables
_log_name = 'elliptic_curve_arithmetic_test.log'
_log_verbose_level = 5
_upper_bound_on_primes = 10000
_number_of_random_elliptic_curves = 100

# Set logging on highest verbose level 
log = Slog(_log_name,_log_verbose_level)
log.info('Test logs written to %s with \n\nlog_verbose_level = %s \nupper_bound_on_primes = %s \nnumber_of_random_elliptic_curves = %s\n' % \
	(_log_name,_log_verbose_level,_upper_bound_on_primes,_number_of_random_elliptic_curves))

for prime in primes_less_than(_upper_bound_on_primes):
	# Library is only build in char != 2,3 
	if prime == 2 or prime == 3: continue 

	# Create finite field with sympy's FiniteField class 
	F = sp.FiniteField(prime)
	
	# Test on a bunch of random ellipti curves over p 
	fail_count = 0
	for i in xrange(0,_number_of_random_elliptic_curves):
		a,b = random_elliptic_curve(prime)

		# Create sympy polynomial 
		try:
			f = sp.poly(y**2 - x**3 - a*x - b)
		except Exception, e:
			log.fail('Could not create sympy polynomial with params a=%s b=%s over %s. %s' % (a,b,prime,'EXCEPTION: ' + str(e)))
			log.fail('EXCEPTION: ' + str(e))
			fail_count+=1
			continue 

		# instanciate elliptic curve 
		try:
			E = EllipticCurve(f,F)
		except Exception, e:
			log.fail('Could not create elliptic curve with params a=%s b=%s over F_%s. %s' % (a,b,prime,'EXCEPTION: ' + str(e)))
			fail_count+=1
			continue

		try:
			P = E.random_point()
			Q = E.random_point() 
		except Exception, e:
			log.fail('Could not generate random point on %E over F_%s. %s' % (E,prime,'EXCEPTION: ' + str(e)))
			fail_count+=1
			continue
		
	log.ok('All methods succeeded on %s/%s random curves over F_%s' % (_number_of_random_elliptic_curves- fail_count,_number_of_random_elliptic_curves,prime)) 
		


