import sys, random

from RandomlyGeneratedTwoGraphForInsplitting import RandomlyGeneratedTwoGraphForInsplitting
from InsplitTwoGraph import InsplitTwoGraph
from RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph
from calculate_h1 import *
from sage.all import matrix, zz, ChainComplex


def calc_homology(g, insplit=False):
  if insplit:
    print("For insplit of random graph")
  else:
    print("For graph")

  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1)
  H1 = C.homology(1)


  res_gemini = calculate_h1_gemini(g.d_1.matrix, g.d_2.matrix)
  H1_gemini=res_gemini['result']
  H1_claude = calculate_h1_claude(g.d_1.matrix, g.d_2.matrix)
  H2_str=g.d_2.ker_str

  return {'H1 sage':H1, 'H1 gemini':H1_gemini, 'H1 claude':H1_claude, 'H2':H2_str}


def calc_homology_and_insplit_homology(n, z):
  print("Generating random graphs...")
  g = RandomlyGeneratedTwoGraphForInsplitting(n,z)

  print("Randomly generated graph full adjacency matrix:")


  print("######## INSPLITTING ###########")
  if g.is_legit:
    gi = InsplitTwoGraph(g, g.v, g.E1, g.E2)

    print("Calculating homology of random graph...")
    g_homology = calc_homology(g)
    print("Calculating homology of the insplit of the random graph...")

    gi_homology= calc_homology(gi, insplit=True)
    print(f"G homology: \nH1 gemini: {g_homology['H1 gemini']}  \nH1 claude: {g_homology['H1 claude']} \nSage H1 {g_homology['H1 sage']} \nH2 {g_homology['H2']} ")
    print(f"In homology: \nH1 gemini: {gi_homology['H1 gemini']}  \nH1 claude: {gi_homology['H1 claude']}  \nSage H1 {gi_homology['H1 sage']} \nH2 {gi_homology['H2']}")

    if (g_homology['H1 gemini'] != gi_homology['H1 gemini']) and (g_homology['H1 claude'] != gi_homology['H1 claude']) and (g_homology['H2'] != gi_homology['H2']):
      fname = './graphs/'
      fname += g_homology['H1 gemini']
      fname += ' ' + gi_homology['H1 gemini']
      fname += ' ' + gi_homology['H2']
      name += 'n_vertices=' + str(g.n)
      name += 'z=' + str(g.z)
      name += '.txt'
      print(f"Diff homologies! Saving to {fname}")

      g.save_to_file('Original graph ' + fname)
      gi.save_to_file('Insplit graph ' + fname)
      exit()


if __name__ == "__main__":
  print()
  if len(sys.argv) < 3:
    print("Usage: python two_graph.py <n: number vertices> <z: roughly upper bounds on number edges>")
    sys.exit(1)

  n = int(sys.argv[1])
  z = int(sys.argv[2])

  calc_homology_and_insplit_homology(n, z)
