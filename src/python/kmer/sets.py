def calc_dictionary_intersection(a, b):
    intersection = {}
    for key in a :
        if key in b :
            intersection[key] = True
    return intersection

def calc_dictionary_symmetric_difference(a, b):
    symmetric_difference = b
    s = a if len(a) < len(b) else b
    for key in a.keys():
        if not key in b:
            symmetric_difference[key] = True
    return symmetric_difference

def print_dictionary_keys(dict):
    print(list(dict.keys()))
