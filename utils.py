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
        if key in types:
            types[key].append(item)
        else:
            types[key] = [item]
    return types.values()

def termsFlattener(termsList):
    flattened = []
    for term in termsList:
        if hasattr(term, 'terms'):
            if len(termsList) == 1 and termsList[0].terms == termsList:
                return termsList
            flattened.extend(termsFlattener(term.terms))
        else:
            flattened.append(term)
    return flattened
