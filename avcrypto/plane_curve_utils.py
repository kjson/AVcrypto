import sympy as sp
import random as rn 

#TODO random_hyper_elliptic_curve

""" Functions for generating elliptic curves for tests """

def random_elliptic_curve(p):
    """
    Creates random values a,b for an elliptic defined over F_p
    Note this generates elliptic curve of the form 
    y**2 - x**3 -x - b == 0, so a,b will be negative. Also 
    both a,b != 0
    """
    a,b = rn.randrange(1,p),rn.randrange(1,p)
    if (-4*a**3 + 27*b**2) % p != 0:
        return a,b 
    else:
        return random_elliptic_curve(p)

def random_hyper_elliptic_curve(p):
    pass