AVcrypto
========

An experimental public key crypto schemes using abelian varieties.

Idea 
------

Implement the arithmetic and functions in selected [algebraic groups](http://en.wikipedia.org/wiki/Algebraic_group) needed to create discrete logarithm based crypto protocals, such as DSA and Diffie-Hellman. 

This was the result of my 4th year honours project for which a report can be found at [ksj.io/AVcrypto.pdf](https://kjson.github.io/AVcrypto.pdf). As a hobby, I'm interested in expanding this library by following the preceeding of this [CRG](https://www.pims.math.ca/scientific/collaborative-research-groups/crg-explicit-methods-abelian-varieties-2015-2018) at PIMS over the next couple years.

For an decent primer on elliptic curve cryptography, check out [http://andrea.corbellini.name/2015/05/17/elliptic-curve-cryptography-a-gentle-introduction/](http://andrea.corbellini.name/2015/05/17/elliptic-curve-cryptography-a-gentle-introduction/).

Algebraic groups so far
------

- Elliptic Curves
- Jacobian of Hyperelliptic curves 

Installation 
------------------
	- git clone https://github.com/kjson/AVcrypto.git 

Dependencies
------------

- Python 2.7 or 3.3+

- [sympy](http://www.sympy.org/en/index.html)

- [numpy](http://www.numpy.org/)






