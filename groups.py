""" Abstract class for groups  """

class Group(object):
    """ Abstract Group """
    def __init__(self, elements,operation):
        self.elements = elements
        self.operation = operation
        self.operationSymbol = '+'

    def __mul__(self,other):
        """ Cartesian product of two groups """
        def operation((x1,y1),(x2,y2)):
            return (self.operation(x1,x2),other.operation(y1,y2))
        return Group((self.elements,other.elements),operation)

    def set_operationSymbol(self,symbol):
        """ Change the symbol for operation """
        self.operationSymbol = symbol

class point(object):
    """ Point in a group,  """
    def __init__(self,point,group):
        assert group.is_elem(point)
        self.group = group
        
    def __add__(self,other):
        """ overload + operator with group operation """
        if self.group.operationSymbol == '+':
            return group.operation(self,other)
        else: 
            Exception("Group operation is %s not +" % self.group.operationSymbol)