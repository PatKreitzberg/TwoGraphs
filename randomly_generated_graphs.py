import sys, random

from InsplitTwoGraph import InsplitTwoGraph
from RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph
from calculate_h1 import *


def calc_homology(g, insplit=False):
  if insplit:
    print("For insplit of random graph")
  else:
    print("For graph")

  res_gemini = calculate_h1_gemini(g.d_1.matrix, g.d_2.matrix)
  H1_gemini=res_gemini['result']
  H1_claude = calculate_h1_claude(g.d_1.matrix, g.d_2.matrix)
  H2_str=g.d_2.ker_str

  return {'H1 gemini':H1_gemini, 'H1 claude':H1_claude, 'H2':H2_str}


def calc_homology_and_insplit_homology(n, z):
  print("Generating random graphs...")
  g = RandomlyGeneratedTwoGraph(n,z)

  print("Randomly generated graph full adjacency matrix:")


  print("######## INSPLITTING ###########")
  if g.is_legit:
    gi = InsplitTwoGraph(g, g.v, g.E1, g.E2)

    print("Calculating homology of random graph...")
    g_homology = calc_homology(g)
    print("Calculating homology of the insplit of the random graph...")
    gi_homology= calc_homology(gi, insplit=True)

    if (g_homology['H1 gemini'] != gi_homology['H1 gemini']) and (g_homology['H1 claude'] != gi_homology['H1 claude']) and (g_homology['H2'] != gi_homology['H2']):
      fname = './graphs/'
      fname += g_homology['H1 gemini']
      fname += ' ' + gi_homology['H1 gemini']
      fname += ' ' + gi_homology['H2']
      name += 'n_vertices=' + str(g.n)
      name += 'z=' + str(g.z)
      name += '.txt'
      g.save_to_file('Original graph ' + fname)
      gi.save_to_file('Insplit graph ' + fname)


if __name__ == "__main__":
  print()
  if len(sys.argv) < 3:
    print("Usage: python two_graph.py <n: number vertices> <z: roughly upper bounds on number edges>")
    sys.exit(1)

  n = int(sys.argv[1])
  z = int(sys.argv[2])

  calc_homology_and_insplit_homology(n, z)
