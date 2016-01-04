import numbers
from functools import reduce
from copy import deepcopy
def isnumeric(n):
    return isinstance(n, numbers.Number)
