from collections import defaultdict
potential_function = defaultdict(tuple)

potential_function[("N", "N")] = 0.2953
potential_function[("N", "P")] = 0.6459
potential_function[("N", "D")] = 0.7114
potential_function[("N", "A")] = 0.6450
potential_function[("N", "H")] = 0.1802
potential_function[("N", "AR")] = 0.0

potential_function[("P", "N")] = 0.6459
potential_function[("P", "P")] = 0.1596
potential_function[("P", "D")] = 0.4781
potential_function[("P", "A")] = 0.7029
potential_function[("P", "H")] = 0.0679
potential_function[("P", "AR")] = 0.1555

potential_function[("D", "N")] = 0.7114
potential_function[("D", "P")] = 0.4781
potential_function[("D", "D")] = 0.5244
potential_function[("D", "A")] = 0.6686
potential_function[("D", "H")] = 0.1453
potential_function[("D", "AR")] = 0.1091

potential_function[("A", "N")] = 0.6450
potential_function[("A", "P")] = 0.7029
potential_function[("A", "D")] = 0.6686
potential_function[("A", "A")] = 0.5478
potential_function[("A", "H")] = 0.2317
potential_function[("A", "AR")] = 0.0770

potential_function[("H", "N")] = 0.1802
potential_function[("H", "P")] = 0.0679
potential_function[("H", "D")] = 0.1453
potential_function[("H", "A")] = 0.2317
potential_function[("H", "H")] = 0.0504
potential_function[("H", "AR")] = 0.0795

potential_function[("AR", "N")] = 0.0
potential_function[("AR", "P")] = 0.1555
potential_function[("AR", "D")] = 0.1091
potential_function[("AR", "A")] = 0.0770
potential_function[("AR", "H")] = 0.0795
potential_function[("AR", "AR")] = 0.1943

print(potential_function)

import pickle

with open('contact_potential_function.pkl', 'wb') as f:
    pickle.dump(potential_function, f)

pkl_file = open('contact_potential_function.pkl', 'rb')

data1 = pickle.load(pkl_file)
print(data1)
print(data1['N','N'])