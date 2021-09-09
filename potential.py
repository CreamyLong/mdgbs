from collections import defaultdict
potential_function = defaultdict(tuple)

potential_function[("NegIonizable", "NegIonizable")]=0.2953
potential_function[("NegIonizable", "PosIonizable")]=0.6459
potential_function[("NegIonizable", "Donor")]=0.7114
potential_function[("NegIonizable", "Acceptor")]=0.6450
potential_function[("NegIonizable", "LumpedHydrophobe")]=0.1802
potential_function[("NegIonizable", "Aromatic")]=0.0

potential_function[("PosIonizable", "NegIonizable")]=0.6459
potential_function[("PosIonizable", "PosIonizable")]=0.1596
potential_function[("PosIonizable", "Donor")]=0.4781
potential_function[("PosIonizable", "Acceptor")]=0.7029
potential_function[("PosIonizable", "LumpedHydrophobe")]=0.0679
potential_function[("PosIonizable", "Aromatic")]=0.1555

potential_function[("Donor", "NegIonizable")] = 0.7114
potential_function[("Donor", "PosIonizable")] = 0.4781
potential_function[("Donor", "Donor")] = 0.5244
potential_function[("Donor", "Acceptor")] = 0.6686
potential_function[("Donor", "LumpedHydrophobe")] = 0.1453
potential_function[("Donor", "Aromatic")] = 0.1091

potential_function[("Acceptor", "NegIonizable")]=0.6450
potential_function[("Acceptor", "PosIonizable")]=0.7029
potential_function[("Acceptor", "Donor")]=0.6686
potential_function[("Acceptor", "Acceptor")]=0.5478
potential_function[("Acceptor", "LumpedHydrophobe")]=0.2317
potential_function[("Acceptor", "Aromatic")]=0.0770

potential_function[("LumpedHydrophobe", "NegIonizable")]=0.1802
potential_function[("LumpedHydrophobe", "PosIonizable")]=0.0679
potential_function[("LumpedHydrophobe", "Donor")]=0.1453
potential_function[("LumpedHydrophobe", "Acceptor")]=0.2317
potential_function[("LumpedHydrophobe", "LumpedHydrophobe")]=0.0504
potential_function[("LumpedHydrophobe", "Aromatic")]=0.0795

potential_function[("Aromatic", "NegIonizable")]=0.0
potential_function[("Aromatic", "PosIonizable")]=0.1555
potential_function[("Aromatic", "Donor")]=0.1091
potential_function[("Aromatic", "Acceptor")]=0.0770
potential_function[("Aromatic", "LumpedHydrophobe")]=0.0795
potential_function[("Aromatic", "Aromatic")]=0.1943

print(potential_function)

import pickle

with open('contact_potential_function.pkl', 'wb') as f:
    pickle.dump(potential_function, f)

pkl_file = open('contact_potential_function.pkl', 'rb')

data1 = pickle.load(pkl_file)
print(data1)