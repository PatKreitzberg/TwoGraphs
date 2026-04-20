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

  def build_matrix(self, graph):
    nv = len(graph.vertices)
    ne = len(graph.edges)
    matrix = [[0]*ne for _ in range(nv)] # matrix of all zeros
    print(f"Matrix of size {ne} by {nv}")

    print("Edge label to index:", [(k,graph.edge_to_index[v.edge_label]) for k,v in graph.edge_label_to_edge.items()])
    print("Vertx label to index:", [(v,graph.vertex_to_index[v]) for v in graph.vertices])

    i = 1       # i = 1 for edges always
    for e_lbl in graph.edges:
      e = graph.edge_label_to_edge[e_lbl]
      edge_index = graph.edge_to_index[e.edge_label]

      ell = 0
      vertex_index = graph.vertex_to_index[e.F(i,ell)]
      matrix[vertex_index][edge_index] -= 1 #-1 because (-1)^(i+ell) = -1
      print(f"F_1^{ell}({e_lbl}) = {e.F(i,ell)}")

      ell = 1
      vertex_index = graph.vertex_to_index[e.F(i,ell)]
      matrix[vertex_index][edge_index] += 1 #+1 because (-1)^(i+ell) = 1
      print(f"F_1^{ell}({e_lbl}) = {e.F(i,ell)}")

    return matrix
