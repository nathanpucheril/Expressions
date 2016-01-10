from Expressions import *
prod = lambda x: reduce(lambda x, y: x*y, x)

def interpolate(*points):
    x, points = X(), set(points)
    return sum([prod([(x - x2)/(x1 - x2) for x2, y2 in points
        if (x1,y1) != (x2,y2)]) * y1 for x1, y1 in points]).simplify()

i = interpolate((-1,1),(0,0),(1,1))
print(i)
