import numbers
from functools import reduce
from copy import deepcopy

isnumeric = lambda n: isinstance(n, numbers.Number)
add = lambda x, y: x + y
reduce = reduce

def list_to_listOfTypes(lst, cmp_fn):
    types = {}
    for item in lst:
        key = cmp_fn(item)
        if types.has_key(key):
            types[key].append(item)
        else:
            types[key] = [item]
    return types.values()
