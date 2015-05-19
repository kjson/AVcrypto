AVcrypto
========

Python crypto library for public key schemes using abelian varieties.

NOTE: To help accommodate the potential use of new groups, I've severly refactored this project. 

Idea 
------

Implement the arithmetic and functions in selected [algebraic groups](http://en.wikipedia.org/wiki/Algebraic_group) needed to create discrete logarithm based crypto protocals, such as DSA and DHE. 

This was the result of my 4th year honours project for which the report can be found [ksj.io/AVcrypto.pdf](http://ksj.io/AVcrypto.pdf). As a hobby, I'm interested in hopefully expanding this library by following the preceeding of this [CRG](https://www.pims.math.ca/scientific/collaborative-research-groups/crg-explicit-methods-abelian-varieties-2015-2018) at PIMS over the next couple years.

For an excellent primer on elliptic curve cryptography, check out [http://andrea.corbellini.name/2015/05/17/elliptic-curve-cryptography-a-gentle-introduction/](http://andrea.corbellini.name/2015/05/17/elliptic-curve-cryptography-a-gentle-introduction/).

Groups
------

- Finite fields of prime order 
- Elliptic Curves
- Jacobian of Hyperelliptic curves 

Long list of todos
------------------

- Add hyperelliptic point counting algorithm 
- Finish SEA algorithm 
- Key exchange on generic groups
- Digital signatures on generi groups 
- Sphinx! 

Dependencies
------------

- Python 2.7 or 3.3+

- [sympy](http://www.sympy.org/en/index.html)

- [numpy](http://www.numpy.org/)






