import sympy as S
from groups import Group   


""" Class for abelian varieties """

class AbelianVariety(Group):
    """
    A complete algebraic abelian group defined as the zero locus of
    a ideal of polynomials over a finite field of prime order.
    Let A r
    def fjsdjfgjgdfgjd 
    epresent an Abelian Variety in all the doc-strings.
    """

    def __init__(self,equations,field,operation):
        
        if field.characteristic() <= 3: # Characteristic != 2,3
            raise NotImplementedError("Only finite fields of prime order are supported.")
  
        symbols = set()
        polynomials = list() # Force sympy Poly class for equations 
        for f in equations:
            symbols = set(symbols | f.atoms(S.Symbol))
            polynomials.append(S.poly(f, symbols, domain=field))

        self.equations = equations
        self.field = field
        self.operation = operation
        self.polynomials = polynomials
        self.symbols = symbols

        # Temp: waiting for solution here: http://stackoverflow.com/questions/29980721/finitefield-object-has-no-attribute-is-commutative
        # Group.__init__(self,S.solve_poly_system(self.polynomials, self.field),operation)
        Group.__init__(self,[],operation)

    def points(self):
        """ Finds all points which satisfy eqautions of A """
        return S.solve_poly_system(self.polynomials, self.field)

    def is_point(self, point):
        """ Verifies that given point belongs to A """

        # Check that point has right number of coordinates
        if len(point) != len(self.symbols):
            raise Exception("Point must be of proper tuple size")

        # f(p) = 0 for all f in I(A)
        for f in self.polynomials:
            if f.eval(point) != 0:
                return False
        return True

    def order(self):
        """ Returns number of points in A """

    def irreducible(self):
        """ Determines if A is irreducible """

    def add(self, p1, p2):
        """ Adds two points on A using defined operation """
        return self.operation(p1, p2)

