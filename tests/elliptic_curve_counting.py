"""
Tests for the arithmetic methods in the EllipticCurve class. Methods are 
tested over a range of finite fields on random elliptic curves. 

methods tested are; order()

"""
from algebraic_groups import EllipticCurve
from slog import Slog
from nt_utils import primes_less_than, random_elliptic_curve

import os
import time

import sympy as sp
import random as rn

# Create sympy variables
x, y = sp.symbols('x, y')

# Set test variables
_log_name = 'elliptic_curve_counting_test.log'
_log_verbose_level = 5
_upper_bound_on_primes = 10000
_number_of_random_elliptic_curves = 100

# Set logging on highest verbose level 
log = Slog(_log_name,_log_verbose_level)
log.info('Test logs written to %s with \n\nlog_verbose_level = %s \nupper_bound_on_primes = %s \nnumber_of_random_elliptic_curves = %s\n' % \
	(_log_name,_log_verbose_level,_upper_bound_on_primes,_number_of_random_elliptic_curves))

test_completion_results = []
for prime in primes_less_than(_upper_bound_on_primes):
	# Library is only build in char != 2,3 
	if prime == 2 or prime == 3: continue 

	# Create finite field with sympy's FiniteField class 
	F = sp.FiniteField(prime)
	
	# Test on a bunch of random ellipti curves over p 
	fail_count = 0
	if _number_of_random_elliptic_curves > prime**2:
		num_curves = prime**2
	else:
		num_curves = _number_of_random_elliptic_curves

	for i in xrange(0,num_curves):
		a,b = random_elliptic_curve(prime)
		f = sp.poly(y**2 - x**3 - a*x - b)
		E = EllipticCurve(f,F)

		"""These method are related to point pointing """

		# Find random points P,Q on E
		try:
			_ = E.order()
		except Exception, e:
			log.fail('COUNTING POINTS: Could not count points for \
			y**2 - x**3 - %sx + %s over F_%s. %s' % (a,b,prime,'EXCEPTION: ' + str(e)))
			fail_count+=1
			continue

	test_completion_results.append((num_curves-fail_count)/num_curves)
	log.info('All methods succeeded on %s/%s random curves over F_%s' % (num_curves- fail_count,num_curves,prime))

test_completion_percentage = sum(test_completion_results)/len(test_completion_results)
log.ok('Tests complete with %s succes rate.' % test_completion_percentage)


