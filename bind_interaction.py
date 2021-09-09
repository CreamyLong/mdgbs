from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from strawberryfields.apps import data, plot, sample, clique

import os
import numpy as np
import networkx as nx
import plotly
import random
from itertools import combinations, permutations
from matrix_utils import check_point_vaild, make_Abbrev_list, make_complete_weight_matrix
import itertools

fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)


def remove_water(m):
  from rdkit.Chem.SaltRemover import SaltRemover
  remover = SaltRemover(defnData="[O]")
  return remover.StripMol(m)


"""输入蛋白质与小分子的距离对"""
def check_pair_match(ligand_pair_distance, protein_pair_distance):

  flexibility_constant_tau = 0
  interaction_distance_epsilon = 0

  if abs(ligand_pair_distance - protein_pair_distance) <= flexibility_constant_tau + 2 * interaction_distance_epsilon:
    # print(abs(ligand_pair_distance - protein_pair_distance))
    return 1
  else:
    return 0

def euclidean_distance(a,b):
  from math import sqrt
  return round(sqrt(sum((a - b)**2 for a, b in zip(a, b))), 1)

def get_feats(mol):
  Donor_feat = factory.GetFeaturesForMol(mol, includeOnly="Donor")
  Acceptor_feat = factory.GetFeaturesForMol(mol, includeOnly="Acceptor")
  NegIonizable_feat = factory.GetFeaturesForMol(mol, includeOnly="NegIonizable")
  PosIonizable_feat = factory.GetFeaturesForMol(mol, includeOnly="PosIonizable")
  Aromatic_feat = factory.GetFeaturesForMol(mol, includeOnly="Aromatic")
  LumpedHydrophobe_feat = factory.GetFeaturesForMol(mol, includeOnly="LumpedHydrophobe")

  feats = Donor_feat + Acceptor_feat + NegIonizable_feat + PosIonizable_feat + Aromatic_feat + LumpedHydrophobe_feat
  return feats

def make_weight_matrix(feats):
  init_matrix = np.zeros((len(feats), len(feats)))
  pharmacophore_list = []

  for i in range(len(feats)):
    p1_pos = feats[i].GetPos()
    p1_type = feats[i].GetFamily()
    # print(p1_type)
    pharmacophore_list.append(p1_type)
    for j in range(len(feats)):
      p2_pos = feats[j].GetPos()
      dis = euclidean_distance(p1_pos, p2_pos)
      init_matrix[i][j] = dis
  processed_weight_matrix = init_matrix
  return processed_weight_matrix, pharmacophore_list


def make_binding_interaction_graph(protein_complete_weight_matrix, ligand_complete_weight_matrix):

  link_list = []

  for p_keys, p_values in protein_complete_weight_matrix.items():
    for l_keys, l_values in ligand_complete_weight_matrix.items():
      res1 = check_pair_match(ligand_pair_distance=p_values, protein_pair_distance=l_values)  # 条件1
      if res1:
        # print(res1)

        tuple_2 = p_keys + l_keys
        # print(list(combinations(tuple_2, 2)))
        allow_pair_list = []
        for pair in list(combinations(tuple_2, 2)):
          res2 = check_point_vaild(pair)   # 条件2
          if res2:
            # print(res2)
            allow_pair_list.append(res2)
        if len(allow_pair_list) == 2:
          #  print(allow_pair_list)
          link_list.append(allow_pair_list)
      else:
        continue

  """筛选逻辑"""
  print("link_list", link_list)
  link_uni_list2 = list(itertools.chain.from_iterable(link_list))
  print("link_uni_list2", link_uni_list2)

  link_count_dict = {}
  for i in link_uni_list2:
    count = link_uni_list2.count(i)
    link_count_dict[i] = count

  # top_link_list = sorted(link_count_dict, key=lambda i: i[1]) #add找到度大的节点
  # print("top_link_list", top_link_list)
  # print(" ")
  # print(sorted(link_count_dict.items(), key=lambda i: i[1], reverse=True)) #add找到度大的节点
  # top_link_list = sorted(link_count_dict.items(), key=lambda i: i[1], reverse=True) #add找到度大的节点
  #
  #
  # print(top_link_list[0][0])
  # link_uni_list = [top_link_list[i][0] for i in range(6)]
  # print(link_uni_list)

  """筛选逻辑"""

  link_uni_list = sorted(list(set(list(itertools.chain.from_iterable(link_list)))))

  # link_uni_list = [link_uni_list[3], link_uni_list[4], link_uni_list[11], link_uni_list[42], link_uni_list[97]]
  # link_uni_list = [link_uni_list[33], link_uni_list[57], link_uni_list[70], link_uni_list[68], link_uni_list[104], link_uni_list[80]] #3rsx

  print("link_uni_list", link_uni_list)

  # print("random")
  # print(random.sample(link_uni_list, 6))  # 结果['a', 'd', 'b', 'f', 'c']，每次运行结果不同。
  # link_uni_list = random.sample(link_uni_list, 5)
  binding_interaction_graph = np.zeros((len(link_uni_list), len(link_uni_list)))
  # binding_interaction_graph = np.zeros(6, 6)

  print("最终图一共有%s个节点" % len(link_uni_list))
  coordinate = {}

  for i in range(len(link_uni_list)):
    for j in range(len(link_uni_list)):

      temp_list = []
      temp_list.append(link_uni_list[i])
      temp_list.append(link_uni_list[j])

      temp_list2 = []
      temp_list2.append(link_uni_list[i])
      temp_list2.append(link_uni_list[j])

      if temp_list in link_list or temp_list2 in link_list:
        binding_interaction_graph[i][j] = 1
        binding_interaction_graph[j][i] = 1

        coordinate[(i, j)] = [link_uni_list[i], link_uni_list[j]]
        coordinate[(j, i)] = [link_uni_list[j], link_uni_list[i]]

        # print(link_uni_list[i], link_uni_list[j])

  print(np.where(binding_interaction_graph == 1))

  print("coordinate", coordinate)

  return binding_interaction_graph, link_uni_list, coordinate


if __name__ == '__main__':
  m = Chem.SDMolSupplier("data/3rsx_ligand.sdf")[0]
  m = Chem.AddHs(m)
  AllChem.EmbedMolecule(m, randomSeed=0xf00d)  # 力场优化

  receptor = Chem.MolFromPDBFile("data/3rsx_pocket.pdb", sanitize=False, removeHs=False)
  receptor = remove_water(receptor)

  l_feats = get_feats(m)  #l represents ligand
  r_feats = get_feats(receptor)   #r represents receptor

  print("ligand's pharmacophore %s points" % len(l_feats))
  print("ligand's pharmacophore %s points name: ", [i.GetFamily() for i in l_feats])

  print("receptor's pharmacophore %s points" % len(r_feats))
  # for i in r_feats:
  #   print(i.GetFamily())

  processed_l_weight_matrix, l_pharmacophore_list = make_weight_matrix(l_feats)
  processed_r_weight_matrix, r_pharmacophore_list = make_weight_matrix(r_feats)


  l_pharmacophore_Abbrev_list = make_Abbrev_list(l_pharmacophore_list, moltype="ligand")
  r_pharmacophore_Abbrev_list = make_Abbrev_list(r_pharmacophore_list, moltype="receptor")

  ligand_complete_weight_matrix = make_complete_weight_matrix(l_pharmacophore_Abbrev_list, processed_l_weight_matrix)
  receptor_complete_weight_matrix = make_complete_weight_matrix(r_pharmacophore_Abbrev_list, processed_r_weight_matrix)

  # print(ligand_complete_weight_matrix, receptor_complete_weight_matrix)
  binding_interaction_graph, link_uni_list, coordinate = make_binding_interaction_graph(ligand_complete_weight_matrix, receptor_complete_weight_matrix) #最终的图， 横轴每个坐标名称

  print(binding_interaction_graph.shape, len(link_uni_list))
  print("link_uni_list", link_uni_list)


  print("binding_interaction_graph", binding_interaction_graph)
  print("coordinate", coordinate)

  np.save("result/3rsx.npy", binding_interaction_graph)

  # TA_graph = nx.Graph(binding_interaction_graph)
  # # print(plot.graph(TA_graph))

  BI = binding_interaction_graph
  BI_graph = nx.Graph(binding_interaction_graph)
  postselected = sample.postselect(list(BI), 2, 10)  # Post-select samples
  samples = sample.to_subgraphs(postselected, BI_graph)  # Convert samples into subgraphs
  shrunk = [clique.shrink(s, BI_graph) for s in samples]  # Shrink subgraphs to cliques
  searched = [clique.search(s, BI_graph, 10) for s in shrunk]  # Perform local search
  clique_sizes = [len(s) for s in searched]
  print("searched ", searched)
  largest_clique = searched[np.argmax(clique_sizes)]  # Identify largest clique found
  print("Largest clique found is = ", largest_clique)

  print("binding sites are")
  for i in list(combinations(largest_clique, 2)):
    print(coordinate[i])
