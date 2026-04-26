import sys, random

from InsplitTwoGraph import InsplitTwoGraph
from RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph
from calculate_h1 import *

def print_adj_matrix(g):
  A = [[0 for _ in range(len(g.vertices))] for _ in range(len(g.vertices))]
  for edge in g.edges:
    A[edge.s][edge.r] += 1
  out = '['
  for r in A:
    for c in r:
      out += str(c) + ' '
    out += ']\n['
  print('\n',out[:-1])


def calc_homology(g, insplit=False):
  if insplit:
    print("For insplit of random graph")
  else:
    print("For graph")

  res_gemini = calculate_h1_gemini(g.d_1.matrix, g.d_2.matrix)
  H1_gemini=res_gemini['result']
  H1_claude = calculate_h1_claude(g.d_1.matrix, g.d_2.matrix)
  H2_str=g.d_2.ker_str # str

  print("H1:")
  print(f"Gemini: H1 = {H1_gemini}")
  print(f"Claude: {H1_claude}")
  print('H2 = ' + H2_str)
  print()


def calc_homology_and_insplit_homology(n, z):
  print("Generating random graphs...")
  g = RandomlyGeneratedTwoGraph(n,z)

  print("Randomly generated graph full adjacency matrix:")
  print_adj_matrix(g)

  print("Calculating homology of random graph...")
  calc_homology(g)

  print("######## INSPLITTING ###########")

  if g.is_legit:
    gi = InsplitTwoGraph(g, g.v, g.E1, g.E2)
    calc_homology(gi, insplit=True)

if __name__ == "__main__":
  print()
  if len(sys.argv) < 3:
    print("Usage: python two_graph.py <n: number vertices> <z: roughly upper bounds on number edges>")
    sys.exit(1)

  n = int(sys.argv[1])
  z = int(sys.argv[2])

  calc_homology_and_insplit_homology(n, z)
