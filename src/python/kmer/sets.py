def calc_dictionary_intersection(a, b):
    intersection = {}
    small = a if len(a.keys()) < len(b.keys()) else b
    large = b if len(a.keys()) < len(b.keys()) else a
    for key in small :
        if key in large :
            intersection[key] = True
    return intersection

def calc_dictionary_union(a, b):
    small = a if len(a.keys()) < len(b.keys()) else b
    large = b if len(a.keys()) < len(b.keys()) else a
    union = large
    for key in small:
        if not key in large:
            union[key] = True
    return union

def calc_dictionary_difference(a, b):
    difference = {}
    for key in a:
        if not key in b:
            difference[key] = True
    return difference

def calc_dictionary_symmetric_difference(a, b):
    symmetric_difference = b
    s = a if len(a) < len(b) else b
    for key in a.keys():
        if not key in b:
            symmetric_difference[key] = True
    return symmetric_differencels

def print_dictionary_keys(dict):
    print(list(dict.keys()))