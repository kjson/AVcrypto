from sympy import abc, poly

""" Modular polynomials for primes up to 23 """

# TODO:
# A Count points 
# A verify irreducibility 
# E count points 
# E Factor 
# H Operation 
# H Count points 
# H Factor 

F,J = symbols('F,J')
x, y = S.symbols('x,y')
f = y ** 2 - x ** 3 - x - 1
F = S.FiniteField(101)
E = EllipticCurve([f], F)
E.schoof()

