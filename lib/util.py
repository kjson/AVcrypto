import math 

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

def easy_factor(n):   
    """ Simple integer factorization """
    return set(x for tup in ([i, n//i] for i in range(1, int(n**0.5)+1) if n % i == 0) for x in tup)
    
def chinese_remainder_theorem(mods,values,len_mods):
    """ Solves modular conguences """
    p = i = prod = 1; sm = 0
    for i in range(len_mods): prod *= n[i]
    for i in range(len_mods):
        p = prod / n[i]
        sm += a[i] * mul_inv(p, n[i]) * p
    return sm % prod





    