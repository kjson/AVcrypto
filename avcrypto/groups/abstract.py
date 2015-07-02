""" Classes for abelian finite groups  """

class Group(object):
    """ Abstract abelian finite group """
    def __init__(self, elements,operation):
        self.elements = points
        self.operation = operation
        
    def __mul__(self,other):
        """ Cartesian product of two groups. Operation is the 
        canonically cartesian operation """
        def operation((x1,y1),(x2,y2)):
            return (self.operation(x1,x2),other.operation(y1,y2))
        return Group((self.points,other.points),operation)

    def is_point(self,point):
        """ Checks if point is in group """   
        if point in self.points:
            return True
        else:
            return False

    def scalar_multiplication(self,scalar,point):
        """ Trivial scalar multiplication of points in group.
        Hopefully a child class of Group has a better implentation. """
        for i in xrange(0,scalar):
            point = point + point
        return point


class Point(object):
    """ Point in an abelian finite group """
    def __init__(self,point,group):
        assert group.is_point(point)
        self.group = group
        self.point = point
        
    def __add__(self,other):
        """ overload + operator with group operation """
        return self.group.operation(self,other)
  
    def __mul__(scalar,self):
        """ Scalar multiplication of a point """  
        return self.group.scalar_multiplication(scalar,self)    

