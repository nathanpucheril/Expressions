import numbers
from functools import reduce
from copy import deepcopy

isnumeric = lambda n: isinstance(n, numbers.Number)
add = lambda x, y: x + y
