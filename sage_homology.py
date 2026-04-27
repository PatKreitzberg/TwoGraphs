import sys, random
import numpy as np

from sage.all import matrix, ZZ, ChainComplex
from RandomlyGeneratedTwoGraphForInsplitting import RandomlyGeneratedTwoGraphForInsplitting
from InsplitTwoGraph import InsplitTwoGraph
from RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph
from TwoGraph import TwoGraph

def calc_homology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1)
  H0 = C.homology(0)
  H1 = C.homology(1)
  H2 = C.homology(2)
  return {'H0':H0, 'H1':H1, 'H2':H2}

def save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2):
    fname = './graphs/'
    fname += 'n_vertices=' + str(g.n)
    fname += ' GH1:' +str(g_H1)
    fname += ' GH2:' +str(g_H2)
    fname += ' GIH1:' +str(gi_H1)
    fname += ' GIH2' +str(gi_H2)

    fname += ' z=' + str(g.z)
    end = ' ' + str(np.random.randint(0,1000000)) + '.txt'
    og_name = fname +' Original graph ' +  end
    ig_name = fname + ' insplit graph ' + end

    if (g_H1 !=  gi_H1) and (g_H2 !=  gi_H2):
        print(f"Diff homologies! Saving to \n{og_name} \nand\n{ig_name}")
        g.save_to_file(og_name)
        gi.save_to_file(ig_name)

        # Append notes
        with open(og_name,'a') as f:
          f.write('#notes\n')
          f.write(' GH1:' +str(g_H1) + "\n")
          f.write(' GH2:' +str(g_H2) + "\n")
          f.write(' GIH1:' +str(gi_H1) + "\n")
          f.write(' GIH2' +str(gi_H2) + "\n")
          f.close()
        exit()

    elif g_H2 !=  gi_H2:
      print("Diff H2!")
    else:
      print("Same homology :(")



def gen_random_graph_and_calc_homology_and_insplit_homology(n, z, runs):
  for run in range(runs):
    print("Run", run)
    g = RandomlyGeneratedTwoGraphForInsplitting(n,z)

    if not g.is_legit:
      print("No insplit vertex found!")
      return

    gi = InsplitTwoGraph(g, g.v, g.E1, g.E2)

    g_homology = calc_homology(g)
    gi_homology = calc_homology(gi)
    g_H0,  g_H1,  g_H2 =  g_homology['H0'],  g_homology['H1'], g_homology['H2']
    gi_H0, gi_H1, gi_H2 = gi_homology['H0'], gi_homology['H1'], gi_homology['H2']

    print(f"G homology: \nH0 {g_H0}, \nH1 {g_H1} \nH2 {g_H2} ")
    print(f"In homology:\nH0 {gi_H0}, \nH1 {gi_H1} \nH2 {gi_H2}")
    save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2)


def premade_graphs(self, og, ig):
  g = TwoGraph(og)
  gi = TwoGraph(ig)



if __name__ == "__main__":
  # Generate random
  if '-random' in sys.argv:
    runs = 1
    if len(sys.argv) == 5:
      runs = int(sys.argv[4])

    n = int(sys.argv[2])
    z = int(sys.argv[3])
    gen_random_graph_and_calc_homology_and_insplit_homology(n, z, runs)
  elif '-files' in sys.argv:
    premade_graphs(sys.argv[2], sys.argv[3])
  else:
    print("Usage:  -random <n: number vertices> <z: roughly upper bounds on number edges> [number of runs]")
    print("Usage: -files <file name of original graph> <file name of insplit graph>")
    sys.exit(1)



'''
Prompt for LaTex:

I am writing a python script which calculates the homology of a
2-graph. I give you the classes for Edge and Commuting
squares. Commuting squares are of the form <range edge> <source edge>
~ < range edge> < source edge>

I have a TwoGraph class which has attributes vertices, edges, and commuting_squares which are lists of those items

if g is an instanec of TwoGraph, I calculate the homology using Sage:
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1)
  H0 = C.homology(0)
  H1 = C.homology(1)
  H2 = C.homology(2)


Your task: Write a python function which uses these to produce a latex
file that writes out the partial matrices and


class CommutingSquare:
  def __init__(self, r1, s1, r2, s2):
    # commuting squares need to have same source and range
    assert r1.r == r2.r
    assert s1.s == s2.s

    self.__setattr__('s1', s1)
    self.__setattr__('r1', r1)
    self.__setattr__('s2', s2)
    self.__setattr__('r2', r2)
    self.__setattr__('path1', (r1, s1))
    self.__setattr__('path2', (r2, s2))
    self.__setattr__('label', '(' + r1.label + ' ' + s1.label + ' ~ ' + r2.label + ' ' + s2.label + ')')

    assert self.path1[0].degree == self.path2[1].degree
    assert self.path1[1].degree == self.path2[0].degree

    self.degree_indices = set([r1.degree, s1.degree, r2.degree, s2.degree])

  def __setattr__(self, name, value):
    if hasattr(self, name):
      print("WARNING: Changed", name)
      raise AttributeError(f"{name} is immutable")
    super().__setattr__(name, value)

  def __eq__(self, other):
    return (self.path1 == other.path1) and (self.path2 == other.path2)

  def __contains__(self, edge):
    assert type(edge) is Edge
    return (edge in self.path1) or (edge in self.path2)

  def __hash__(self):
    return hash(self.label)

  def __str__(self):
    return self.label

  def __getitem__(self, i):
    assert i < 4
    if i <= 1:
      return self.path1[i]
    return self.path2[i-2]


  def F(self, i, ell):
    # i is degree index
    assert i in self.degree_indices
    assert ell in [0,1]

    ell_index = (ell+1)%2 # Because F_i^0 should be the range edge

    if self.path1[ell_index].degree == i:
      return self.path1[ell_index]

    return self.path2[ell_index]


class Edge:
  def __init__(self, label, s, r, degree=None):
    # label is string name of edge
    # s is source vertex
    # r is range vertex

    self.label = label
    self.s = s # source of edge
    self.r = r # range of edge
    self.degree = degree
    self.range_of_commuting_squares = set()
    self.commuting_squares = set()

  def __contains__(self, v):
    return (v == self.s) or (v == self.r)
  def __str__(self):
    return self.label
  def __eq__(self, other):
    return self.label == other.label
  def __hash__(self):
    return hash(self.label)
  def __lt__(self, other):
    return self.label < other.label

  def F(self, i, ell):
    assert i == 1
    assert ell in [0,1]
    if ell == 0:
      return self.r
    if ell == 1:
      return self.s


'''
