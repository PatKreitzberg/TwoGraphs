from sympy import Matrix,pprint
from sympy.matrices.normalforms import smith_normal_form


class PartialOneMatrix:
  def __init__(self, graph):
    self.matrix = self.build_matrix(graph)

  def __str__(self):
    out = ''
    for row in self.matrix:
      out += str(row) + '\n'
    return out

  def increment(self, i, ell):
    if (i + ell)%2 == 0:
      return 1
    return -1

  def build_matrix(self, graph):
    nv = len(graph.vertices)
    ne = len(graph.edges)
    matrix = [[0]*ne for _ in range(nv)] # matrix of all zeros

    i = 1       # i = 1 for edges always
    for e_lbl in graph.edges:
      e = graph.edge_label_to_edge[e_lbl]
      edge_index = graph.edge_to_index[e.edge_label]
      for ell in [0,1]:
        vertex_index = graph.vertex_to_index[e.F(i,ell)]
        matrix[vertex_index][edge_index] += self.increment(i,ell)

    return matrix


class PartialTwoMatrix:
  def __init__(self, graph):
    self.matrix = self.build_matrix(graph)

  def increment(self, i, ell):
    if (i + ell)%2 == 0:
      return 1
    return -1

  def build_matrix(self, graph):
    ne = len(graph.edges)
    ncs = len(graph.commuting_squares)
    matrix = [[0]*ncs for _ in range(ne)]

    for cs in self.commuting_squares:
      cs_index = graph.commuting_square_to_index[cs]
      for i in [1,2]:
        for ell in [0,1]:
          edge_index = graph.edge_to_index[cs.F(i, ell)]
          matrix[edge_index][cs_index] += self.increment(i,ell)
    return matrix
