def calc_dictionary_intersection(a, b):
    intersection = {}
    small = a if len(a.keys()) < len(b.keys()) else b
    large = b if len(a.keys()) < len(b.keys()) else a
    for key in small :
        if key in large :
            intersection[key] = True
    return intersection

def calc_dictionary_symmetric_difference(a, b):
    symmetric_difference = b
    s = a if len(a) < len(b) else b
    for key in a.keys():
        if not key in b:
            symmetric_difference[key] = True
    return symmetric_differencels

def print_dictionary_keys(dict):davisss
    print(list(dict.keys()))