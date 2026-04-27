import sys, random
import numpy as np

from RandomlyGeneratedTwoGraphForInsplitting import RandomlyGeneratedTwoGraphForInsplitting
from InsplitTwoGraph import InsplitTwoGraph
from RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph
from TwoGraph import TwoGraph

def print_graph(g):
  g_d1_str = g.d_1.latex()
  g_d2_str = g.d_2.latex()
  g.d_1.calc_img()
  g.d_1.calc_ker()

  print(g_d1_str)
  print()
  print(g_d2_str)
  print()

  print("Img", g.d_1.img_str_items)
  print()
  print("Ker", g.d_1.ker_str_items)

def premade_graphs(og, ig):
  g = TwoGraph(og)
  gi = TwoGraph(ig)

  print_graph(g)
  print_graph(gi)


if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: <file name of original graph> <file name of insplit graph>")
    sys.exit(1)

premade_graphs(sys.argv[1], sys.argv[2])
