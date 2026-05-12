import sys,re
from abc import ABC, abstractmethod

import numpy as np

from python.Edge import Edge
from python.CommutingSquare import CommutingSquare

from sage.all import vector, matrix, ZZ


class Graph(ABC):
  def __init__(self, k):
    self.k = k
    self.edges, self.commuting_squares  = self.build(k)
    self.k_to_cubes = dict()
    self.k_to_cubes[k] = self.generate_kcubes()
    self.finish_kcubes()
    self.k_to_cubes[k-1] = self.k_minus_1_cubes()
    self.build_dk()
    self.print_dk()
    #self.print_kcubes()

  @abstractmethod
  def build(self, k):
    pass

  @abstractmethod
  def generate_kcubes(self):
    pass

  def strip(self, edge):
    # Removes the _i from an edge
    return re.sub(r'_\d+', '', edge)

  def add_edge_index(self, edge, index):
    return edge[0] + '_' + str(index) + edge[1:]

  def print_kcubes(self):
    print("K cubes for", self.label, 'with k=', self.k)
    for k,cubes in self.k_to_cubes.items():
      print(f"k={k}")
      for cube in cubes:
        print(cube)
      print()

  def finish_kcube(self, cube):
    for j, elt in enumerate(cube):
      elt = elt[0] + '_' + str(j+1) + elt[1:]
      cube[j] = elt
    return tuple(cube)

  def finish_kcubes(self):
    for k,cubes in self.k_to_cubes.items():
      for i in range(len(cubes)):
        cubes[i] = self.finish_kcube(list(cubes[i]))

  def move_j_index_to_end(self, ell, j, cube):
    '''
    Need to rewrite cubes to move the jth index to the end. May be the first spot, may be the last

    j starts at 0

    ell=0 means move to first spot
    ell=1 means move to last spot
    '''
    newcube = []
    edge_indices = list(range(len(cube)))

    if ell == 0:
      ei = [j] + edge_indices[:j] + edge_indices[(j+1):]
    else:
      ei = edge_indices[:j] + edge_indices[(j+1):] + [j]

    pattern = r'_\d+(?=\^|$)'
    for index,elt in zip(ei, cube):
      newcube.append(re.sub(pattern, f'_{str(index+1)}', elt))

    newcube = tuple(newcube)
    return newcube

  def F(self, ell, j, cube):
    # Get face from kcube
    cube = self.move_j_index_to_end(ell, j, cube)

    assert ell in [0,1]
    if ell == 0:
      return cube[1:]
    return cube[:-1]

  def k_minus_1_cubes(self):
    cubes = set()
    # only build necessary kcubes
    for cube in self.k_to_cubes[self.k]:
      for j in range(self.k):
        for ell in [0,1]:
          cubes.add(self.F(ell, j, cube))
    return list(cubes)

  def ker_rank(self):
    M = matrix(ZZ, self.dk)
    return M.right_kernel().rank()

  def print_dk(self):
    print("\ndk matrix")
    for r in self.dk:
      print(r)
    print()

  def build_dk(self):
    '''Boundary map'''

    # dk has one column per kcube
    #            one row      per (k-1)cube

    self.dk = [[0]*len(self.k_to_cubes[self.k]) for _ in range(len(self.k_to_cubes[self.k-1]))]

    k_cube_to_index               = {cube:i for i, cube in enumerate(self.k_to_cubes[self.k])}
    k_minus_1_cube_to_index = {cube:i for i, cube in enumerate(self.k_to_cubes[self.k-1])}

    for k_cube, k_cube_index in k_cube_to_index.items():
      for j in range(self.k):
        for ell in [0,1]:
          sign = 1
          if (ell+j)%2 == 1:
            sign = -1

          k_minus_1_cube = self.F(ell, j, k_cube)
          k_minus_1_cube_index = k_minus_1_cube_to_index[k_minus_1_cube]
          row = k_minus_1_cube_index
          col  = k_cube_index
          self.dk[row][col] += sign


class SingleVertexKGraph(Graph):
  '''A single vertex k-graph will contain k-cubes which can be
  represented as tuples (i_1, i_2,...) with all i_j in {1,2}. This is
  because the commuting squares will be:
  ei^1 ej^1 ~ ej^1 ei^1
  ei^1 ej^2 ~ ej^1 ei^2
  ei^2 ej^1 ~ ej^2 ei^1
  ei^2 ej^2 ~ ej^2 ei^2

  where ei are degree i and ej are degree j. In each
  commuting square, the superscript stays the same while the edge
  degree designation shifts. So, a kcube is determined by
  only the index: If k=4 then
  (1,1,2,1) is  e1^1 e2^1 e3^2 e4^1
  (2,2,1,1) is  e1^2 e2^2 e3^1 e4^1
  and so on.
  '''
  def __init__(self, k):
    self.label = '(Single vertex)' + str(np.random.randint(10000))
    super().__init__(k)


  def generate_kcubes(self, kcube=None, kcubes=None):
    if kcube is None:
        kcube = []
    if kcubes is None:
        kcubes = []

    labels = ('a', 'b')
    # tupe is list with either 1 or 2
    if len(kcube) < self.k:
      kcubes = self.generate_kcubes(kcube+[labels[0]], kcubes=kcubes)
      return      self.generate_kcubes(kcube+[labels[1]], kcubes=kcubes)
    else:
      kcubes.append(tuple(kcube))
      return kcubes

  def build(self, k):
    edges,commuting_squares = [],[]
    for k1 in range(1,k+1):
      for k2 in range(k1+1, k+1):
        a_1 = Edge('a_' + str(k1), 'v', 'v', degree=k1)
        b_1 = Edge('b_' + str(k1), 'v', 'v', degree=k1)
        a_2 = Edge('a_' + str(k2), 'v', 'v', degree=k2)
        b_2 = Edge('b_' + str(k2), 'v', 'v', degree=k2)
        edges.append(a_1)
        edges.append(b_1)
        edges.append(a_2)
        edges.append(b_2)
        commuting_squares.append( CommutingSquare(a_1, a_2, a_2, a_1) )
        commuting_squares.append( CommutingSquare(a_1, b_2, a_2, b_1) )
        commuting_squares.append( CommutingSquare(b_1, a_2, b_2, a_1) )
        commuting_squares.append( CommutingSquare(b_1, b_2, b_2, b_1) )
    return edges, commuting_squares


class InsplitSingleVertexGraph(Graph):
  def __init__(self, g):
    self.parentg = g
    self.label = '(Insplit graph)'  + str(np.random.randint(10000))
    super().__init__(g.k)

  def build(self, k):
    self.label_to_edge = {}
    commuting_squares = []
    edges = []
    for k1 in range(1,k+1):
      for k2 in range(k1+1, k+1):
        a1_1 = Edge('a^1_' + str(k1), 'v^1', 'v^1', degree=k1)
        a1_2 = Edge('a^1_' + str(k2), 'v^1', 'v^1', degree=k2)
        a2_1 = Edge('a^2_' + str(k1), 'v^2', 'v^1', degree=k1)
        a2_2 = Edge('a^2_' + str(k2), 'v^2', 'v^1', degree=k2)
        b1_1 = Edge('b^1_' + str(k1), 'v^1', 'v^2', degree=k1)
        b1_2 = Edge('b^1_' + str(k2), 'v^1', 'v^2', degree=k2)
        b2_1 = Edge('b^2_' + str(k1), 'v^2', 'v^2', degree=k1)
        b2_2 = Edge('b^2_' + str(k2), 'v^2', 'v^2', degree=k2)

        self.label_to_edge[a1_1.label] = a1_1
        self.label_to_edge[a1_2.label] = a1_2
        self.label_to_edge[a2_1.label] = a2_1
        self.label_to_edge[a2_2.label] = a2_2
        self.label_to_edge[b1_1.label] = b1_1
        self.label_to_edge[b1_2.label] = b1_2
        self.label_to_edge[b2_1.label] = b2_1
        self.label_to_edge[b2_2.label] = b2_2

        edges.append(a1_1)
        edges.append(b1_1)
        edges.append(a1_2)
        edges.append(b1_2)
        edges.append(a2_1)
        edges.append(b2_1)
        edges.append(a2_2)
        edges.append(b2_2)
        commuting_squares.append(CommutingSquare(a1_1, a2_2, a1_2, a2_1))
        commuting_squares.append(CommutingSquare(a1_1, a1_2, a1_2, a1_1))
        commuting_squares.append(CommutingSquare(a2_1, b2_2, a2_2, b2_1))
        commuting_squares.append(CommutingSquare(a2_1, b1_2, a2_2, b1_1))
        commuting_squares.append(CommutingSquare(b1_1, a2_2, b1_2, a2_1))
        commuting_squares.append(CommutingSquare(b1_1, a1_2, b1_2, a1_1))
        commuting_squares.append(CommutingSquare(b2_1, b1_2, b2_2, b1_1))
        commuting_squares.append(CommutingSquare(b2_1, b2_2, b2_2, b2_1))
    return edges, commuting_squares

  def build_kcube(self, parent, child):
    '''
    parent will be of the from a_1, a_2, b_3, b_4, a_5
    Child will be the first element but insplit so starts with either a_5^1 or a_5^2
    '''
    next_from_prev = {
      'a':{'b^1':'a^2','b^2':'a^2','a^1':'a^1','a^2':'a^1'},
      'b':{'b^1':'b^2','b^2':'b^2','a^1':'b^1','a^2':'b^1'}
    }

    for v in parent[:-1][::-1]: # go backwards and skip first one
      prev = child[-1]
      prev = prev[0] + prev[-2:]
      n = next_from_prev[v[0]][prev]
      assert self.label_to_edge[n + '_1'].s == self.label_to_edge[prev + '_1'].r
      child.append(n)

    return child

  def generate_kcubes(self):
    kcubes = []
    for kcube in self.parentg.k_to_cubes[self.k]:
      # Of the form ('a','b','a',...)

      # For each kcube in g
      # Add two kcubes

      child = self.strip(kcube[-1]) # will look like a_3

      kcubes.append(tuple( self.build_kcube(kcube, [child + '^1'] )[::-1]))
      kcubes.append(tuple( self.build_kcube(kcube, [child+ '^2'] )[::-1]))
    return kcubes

def run(k):
  g = SingleVertexKGraph(k)
  gi = InsplitSingleVertexGraph(g)
  print("Lambda Kernel rank:", g.ker_rank())
  print("Insplit    Kernel rank:", gi.ker_rank())

if __name__ == "__main__":
  run(int(sys.argv[1]))
