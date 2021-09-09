


def find_letter(name:str):
  import re
  pattern = re.compile(r'[a-zA-Z]')   # 查找字母
  result = pattern.findall(name)
  return "".join(result)
find_letter("Ar1")



def check_point_vaild(pair):

    allow_pair = {
        "D": "a",
        "d": "A",
        "a": "D",
        "A": "d",
        "n": "P",
        "N": "p",
        "p": "N",
        "P": "n",
        "H": "h",
        "h": "H",
        "Ar": "ar",
        "ar": "Ar"
    }

    p1 = pair[0]
    p2 = pair[1]
    pl1 = find_letter(p1)
    pl2 = find_letter(p2)
    if pl2 == allow_pair[pl1]:
        return pair
    else:
        return 0

def make_complete_weight_matrix(pharmacophore_Abbrev_list, processed_matrix):
  from collections import defaultdict
  make_complete_weight_matrix = defaultdict(tuple)
  length = processed_matrix.shape[0]
  for i in range(length):
    # if i == 0: continue
    pair_part1 = pharmacophore_Abbrev_list[i]
    for j in range(i+1, length):
      # if j == 0: continue
      pair_part2 = pharmacophore_Abbrev_list[j]
      # print(i,j)
      # print(pair_part1,pair_part2)
      make_complete_weight_matrix[(pair_part1, pair_part2)] = processed_matrix[i][j]
  return make_complete_weight_matrix

def make_Abbrev_list(pharmacophore_list: list, moltype: str = "ligand") -> list:
    Donor = 0
    Acceptor = 0
    NegIonizable = 0
    PosIonizable = 0
    Aromatic = 0
    LumpedHydrophobe = 0

    pharmacophore_Abbrev_list = []

    if moltype == "ligand":
        for i in pharmacophore_list:
            if i == "Donor":
                name = "D" + str(Donor + 1)
                pharmacophore_Abbrev_list.append(name)
                Donor = Donor + 1
            elif i == 'Acceptor':
                name = "A" + str(Acceptor + 1)
                pharmacophore_Abbrev_list.append(name)
                Acceptor = Acceptor + 1
            elif i == 'NegIonizable':
                name = "N" + str(NegIonizable + 1)
                pharmacophore_Abbrev_list.append(name)
                NegIonizable = NegIonizable + 1
            elif i == 'PosIonizable':
                name = "P" + str(PosIonizable + 1)
                pharmacophore_Abbrev_list.append(name)
                PosIonizable = PosIonizable + 1
            elif i == 'Aromatic':
                name = "Ar" + str(Aromatic + 1)
                pharmacophore_Abbrev_list.append(name)
                Aromatic = Aromatic + 1
            elif i == 'LumpedHydrophobe':
                name = "H" + str(LumpedHydrophobe + 1)
                pharmacophore_Abbrev_list.append(name)
                LumpedHydrophobe = LumpedHydrophobe + 1
            else:
                raise Exception("the pharmacophore name is not existed! ")
        return pharmacophore_Abbrev_list

    elif moltype == "receptor":
        for i in pharmacophore_list:
            if i == "Donor":
                name = "d" + str(Donor + 1)
                pharmacophore_Abbrev_list.append(name)
                Donor = Donor + 1
            elif i == 'Acceptor':
                name = "a" + str(Acceptor + 1)
                pharmacophore_Abbrev_list.append(name)
                Acceptor = Acceptor + 1
            elif i == 'NegIonizable':
                name = "n" + str(NegIonizable + 1)
                pharmacophore_Abbrev_list.append(name)
                NegIonizable = NegIonizable + 1
            elif i == 'PosIonizable':
                name = "p" + str(PosIonizable + 1)
                pharmacophore_Abbrev_list.append(name)
                PosIonizable = PosIonizable + 1
            elif i == 'Aromatic':
                name = "ar" + str(Aromatic + 1)
                pharmacophore_Abbrev_list.append(name)
                Aromatic = Aromatic + 1
            elif i == 'LumpedHydrophobe':
                name = "h" + str(LumpedHydrophobe + 1)
                pharmacophore_Abbrev_list.append(name)
                LumpedHydrophobe = LumpedHydrophobe + 1
            else:
                raise Exception("the pharmacophore name is not existed! ")
        return pharmacophore_Abbrev_list