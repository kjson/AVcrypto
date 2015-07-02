import math
import random as rn
from fractions import gcd


def extended_gcd(aa, bb):
    """ Computes the extended euclidean algorithm """
    lastremainder, remainder = abs(aa), abs(bb)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = remainder, \
            divmod(lastremainder, remainder)
        x, lastx = lastx - quotient*x, x
        y, lasty = lasty - quotient*y, y
    return lastremainder, lastx * (-1 if aa < 0 else 1),\
        lasty * (-1 if bb < 0 else 1)


def mod_inv(a, m):
    """ Modular inverse """
    g, x, y = extended_gcd(a, m)
    if g != 1:
		raise ValueError
    return x % m

def chinese_remainder_theorem(n, a, lena):
    p = i = prod = 1; sm = 0
    for i in range(lena): prod *= n[i]
    for i in range(lena):
        p = prod / n[i]
        sm += a[i] * mod_inv(p, n[i]) * p
    return sm % prod


def primes_less_than(n):
    """ Find all primes less than a given bound """
    sieve = [True] * n
    for i in xrange(3,int(n**0.5)+1,2):
        if sieve[i]:
            sieve[i*i::2*i]=[False]*((n-i*i-1)/(2*i)+1)
    return [2] + [i for i in xrange(3,n,2) if sieve[i]]


def is_probable_prime(N):
    """ Miller-Rabin primality test on an integer N """
    # special cases
    if N == 2:
        return True
    elif N==1:
        return False
    elif N % 2 == 0:
        return False

    # write N-1 as 2**s * d
    # repeatedly try to divide N-1 by 2
    s = 0
    d = N-1
    while True:
        quotient, remainder = divmod(d, 2)
        if remainder == 1:
            break
        s += 1
        d = quotient
    assert(2**s * d == N-1)

    _mrpt_num_trials = 20
    # test the base a to see whether it is a witness for the compositeness of N
    def try_composite(a):
        if a**d % N == 1:
            return False
        for i in range(s):
            if a**(2**i*d) % N == N-1:
                return False
        return True # N is definitely composite

    for i in range(_mrpt_num_trials):
        a = rn.randrange(2, N)
        if try_composite(a):
            return False

    return True


def lenstra(N):
    """ Lenstra's elliptic curve factorization method """

    # If N is negative
    if N < 0:
        return [-1] + lenstra(-N)

    # If N is prime
    if 0 <= N <= 3 or is_probable_prime(N):
        return [N]

    # Initialize two random points mod N
    x0,y0 = rn.randrange(1, N),rn.randrange(1, N)

    # Bound for the number of trials
    bound = int(math.sqrt(N))

    for a in xrange(1, N):
        # build elliptic curve using random points
        b = y0**2 - x0**3 - a*x0

        # Check if the curve is singular
        if 4*a**3 - 27*b**2 == 0:
            next

        d = gcd(2*y0, N)
        if d != 1:
            return lenstra(d) + lenstra(N/d)

        # Initial double of the point (x0,y0)
        s = (3*x0**2 + a) * mod_inv(2*y0,N)
        x,y = (s**2 - 2*x0,s*(3*x0 - s**2)-y0)

        # Search for non-trivial gcd's
        for k in xrange(1,bound):
            for i in xrange(1,k):
                d = gcd(x-x0,N)
                if d != 1:
                    return lenstra(d) + lenstra(N/d)
                else:
                    s = (y-y0) * mod_inv(x-x0,N)
                    y = (s*(2*x + x0 - s**2)-y) % N
                    x = (s**2 - x -x0) % N

def compute_modular_polynomial(l):
    """ computes the lth modular polynomial using ISOGENIE VOLCANOS """


def legendre(a,p):
    """ The legendre symbol of a mod p """

    # p must be an odd prime
    if not is_probable_prime(p) or p <=2:
        raise Exception("p must be an odd prime")

    def f(q,p):
        """ Legendre symbol where q is prime """

        # all the base cases
        if q == 1:
            return 1
        if q % p ==0:
            return 0
        if q == -1:
            if p % 4 == 1:
                return 1
            else:
                return -1
        if q == 2:
            if p % 8 == 1 or p % 8 == 7:
                return 1
            if p % 8 == 3 or p % 8 == 5:
                return -1

        # Recursive step
        if q > p:
            return f(q % p,p)
        if q % 4 == 1 or p % 4 == 1:
            return f(p,q)
        else:
            return -f(p,q)

    # Because the legendre symbol is multiplicative
    s = 1
    for x in lenstra(a%p):
        s *= f(x,p)
    return s




