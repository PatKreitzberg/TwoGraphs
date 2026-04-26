import sys

from TwoGraph import TwoGraph
from InsplitTwoGraph import InsplitTwoGraph
from calculate_h1 import calculate_h1_gemini, calculate_h1_claude


def calc_homology(g, insplit=False):
  if insplit:
    print("For insplit of graph:" + path)
  else:
    print("For graph:" + path)

  res_gemini = calculate_h1_gemini(g.d_1.matrix, g.d_2.matrix)
  H1_gemini=res_gemini['result']
  H1_claude = calculate_h1_claude(g.d_1.matrix, g.d_2.matrix)
  H2_str=g.d_2.ker_str # str

  print("H1:")
  print(f"Gemini: H1 = {H1_gemini}")
  print(f"Claude: {H1_claude}")
  print('H2 = ' + H2_str)
  print()


path='graphs/EFGGGP_Figure3.txt'
g = TwoGraph(load_from=path)
calc_homology(g)

gi = InsplitTwoGraph(g)
