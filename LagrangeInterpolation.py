from Expressions import *
prod = lambda x: reduce(lambda x, y: x*y, x)

def interpolate(points):
    """Takes in a list of tuples of x,y cordinates.

    >>> coords = lambda y: [(x, y(x)) for x in range(y.degree+1)]
    >>> interpolate(*coords(P('x^5 + 3x^3 + x')))
    x^5 + 3x^3 + x
    >>> interpolate(*coords(P('7x^4 + 5x^3 + 9x^2 + 3')))
    7x^4 + 5x^3 + 9x^2 + 3
    """
    points = set(points)
    x = X()
    return sum([prod([(x - x2)/(x1 - x2)for x2, y2 in points
        if (x1,y1) != (x2,y2)]) * y1 for x1, y1 in points]).simplify()
