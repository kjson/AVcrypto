import math 
from fractions import gcd
import random 

def extended_gcd(aa, bb):
    """ Computes the extended euclidean algorithm """
    lastremainder, remainder = abs(aa), abs(bb)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(lastremainder, remainder)
        x, lastx = lastx - quotient*x, x
        y, lasty = lasty - quotient*y, y
    return lastremainder, lastx * (-1 if aa < 0 else 1), lasty * (-1 if bb < 0 else 1)

def modInv(a, m):
    """ Modular inverse """
    g, x, y = extended_gcd(a, m)
    if g != 1:
		raise ValueError
    return x % m


def primes_less_than(n):
    """ Find all primes less than a given bound """
    sieve = [True] * n

    for i in xrange(3,int(n**0.5)+1,2):
        if sieve[i]:
            sieve[i*i::2*i]=[False]*((n-i*i-1)/(2*i)+1)

    return [2] + [i for i in xrange(3,n,2) if sieve[i]]


def is_probable_prime(N):
    """ Miller-Rabin primality test on an integer N """
    assert N >= 2
    # special case 2
    if N == 2:
        return True
    # ensure N is odd
    if N % 2 == 0:
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

    # This needs to be set, but 10 corresponds to false possitive probablility = 1/2**10
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
        a = random.randrange(2, N)
        if try_composite(a):
            return False
 
    return True # no base tested showed N as composite

def lenstra_elliptic_curve_factor(N):
    """ Lenstra's elliptic curve factoring method """

    # Base case 
    if N == 1:
        return []

    # Cant factor a prime! 
    if N == 2 or is_probable_prime(N):
        return [N]

    # Initialize two random integers mod N 
    x0, y0 = random.randrange(1, N), random.randrange(1, N)

    # List of factors to be returned 
    factors = list() 

    # 
    bound = int(math.sqrt(N))

    for a in xrange(2,N):
        # Build curve out of random points
        b = y0**2 - x0**3 - a*x0

        # Check curve is not singular 
        if 4*a**3 - 27*b**2 ==0:
            next

        # Initially double point 
        s = (3*x0**2 + a) 
        (x,y) = (s**2 - 2*x0, s*((s**2 - 2*x0) - x0) - y0)

        # Keep adding points until gcd(x-x0,N) != 1
        for k in xrange(2,bound):
            for i in xrange(0,math.factorial(k)):
                d = gcd(x- x0,N)
                if d != 1:
                    return lenstra_elliptic_curve_factor(int(d)) + lenstra_elliptic_curve_factor(int(N/d))
                else:
                    s = (y - y0) * modInv(x - x0,N)
                    x = s**2 - x - x0  
                    y = - y + s * (s**2 - x - x0 - x)
    

def chinese_remainder_theorem(mods,values,len_mods):
    """ Solves modular conguences """
    p = i = prod = 1; sm = 0
    for i in range(len_mods): prod *= n[i]
    for i in range(len_mods):
        p = prod / n[i]
        sm += a[i] * mul_inv(p, n[i]) * p
    return sm % prod


def legendre(a,p):
    """ Computes the legendre sybol of a mod p """

    # p must be an odd prime 
    if not is_probable_prime(p) or p <=2:
        raise Exception("p must be prime")

    # result 
    s = 1

    def f(q,p):
        """ Cqlculqtes legendre symbol where q is prime """

        # qll the bqse cqses 
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
    for x in lenstra_elliptic_curve_factor(a):
        s *= f(x,p)

    return s 





    