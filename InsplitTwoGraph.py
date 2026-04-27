from TwoGraph import TwoGraph
from Edge import Edge
from CommutingSquare import CommutingSquare

class InsplitTwoGraph(TwoGraph):
  def __init__(self, g, v, E1, E2):
    super().__init__()
    # og_graph is the graph we are going to insplit
    # v is the vertex at which we insplit
    self.E1 = E1
    self.E2 = E2
    self.v = v

    assert len(self.E1&self.E2) == 0
    assert len(self.E1) > 0
    assert len(self.E2) > 0

    self.vertices, self.edges, self.commuting_squares = self.insplit(g,v)
    self.n = len(self.vertices)
    self.calculate_boundary_matrices()

  def insplit(self, g, v):
    #Add vertices to graph
    new_vertices = [str(v)+'^1', str(v)+'^2']
    vertices_as_set = self.add_insplit_vertices(v, g.vertices, new_vertices)

    # MODIFY EDGES
    s_inv_e = set(g.source_inverse_of_vertex(v))
    edges, edge_to_children_edges = self.add_insplit_edges(set(g.range_inverse_of_vertex(v)), s_inv_e, set(g.edges), new_vertices, v)

    # COMMUTING SQUARES
    commuting_squares = self.add_insplit_commuting_squares(v, g.commuting_squares, edge_to_children_edges)
    return list(vertices_as_set), edges, commuting_squares

  def add_insplit_source_edges(self, edges, s_inv_e, v, new_vertices):
    # Make function:
    #  e -> [e^1, e^2] if s(e) = v
    #  else e -> e

    edges -= s_inv_e
    edge_to_children_edges= {e:[e] for e in edges} # these edges are not modified at all

    # ADD NEW SOURCE EDGES
    for e in s_inv_e:
      edge_to_children_edges[e] = list()

      for j in [0,1]:
        # range remains the same
        r = e.r

        # make new edges e^1, e^2
        vertex_label = str(j+1)
        label = e.label + "^" + vertex_label
        s = new_vertices[j]

        new_e = Edge(label, s, r, degree=e.degree)

        edges.add(new_e)
        edge_to_children_edges[e].append(new_e)

    return edges, edge_to_children_edges


  def add_insplit_commuting_squares(self, v, old_commuting_squares, edge_to_children_edges):
    # we have af ~ eb if
    # 1. The parent squares in the original graph commute
    # 2. The range(af) = range(eb)
    # 3. The source(af) = source(eb)

    # If an edge was not involved in insplit then

    new_commuting_squares = set()

    for cs in old_commuting_squares:
      (a,e),(f,b) = cs.path1, cs.path2
      for a_ in edge_to_children_edges[a]:
          for e_ in edge_to_children_edges[e]:
              for f_ in edge_to_children_edges[f]:
                  for b_ in edge_to_children_edges[b]:
                    if self.check_commuting_square_is_valid(a_, e_, f_, b_):
                      new_commuting_squares.add(CommutingSquare(a_, e_, f_, b_))

    return list(new_commuting_squares)

  def check_commuting_square_is_valid(self, a, e, f, b):
    # a e ~ f b
    # need:
    # 1. We have two paths
    #     range of e is source of a
    #     range of b is source of f
    # 2. The two paths have same source and range:
    #     source of e is source of b
    #     range of a is range of f
    return (a.s == e.r) and (a.r == f.r) and (e.s == b.s) and (b.r == f.s)


  def add_insplit_vertices(self, v, vertices, new_vertices):
    '''
    Just removes old v and adds v^1 and v^2
    '''
    vertices = set(vertices)
    vertices.remove(v)
    vertices |= set(new_vertices)
    return vertices

  def add_insplit_edges(self, r_inv_e, s_inv_e, edges, new_vertices, v):
    '''
    We have a specific vertex v
    v is replaced with v^1, v^2
    Edges with RANGE V:
      If e has range v then we have to change its range from v to either v^1 or v^2.
      Partition r^{-1}(v) into E1 and E2 such that E1, E2 are nonempty, E1&E2 is empty, and E1 U E2 = r^{-1}(v)
      If e in E1 then range(e) should become V1. Similar if e in E2.

    Edges with SOURCE V:
      If s(e) = v then need to duplicate e and
    '''

    edges = self.add_insplit_range_edges(r_inv_e, new_vertices, edges, v)
    edges, edge_to_children_edges = self.add_insplit_source_edges(edges, s_inv_e, v, new_vertices)
    return edges, edge_to_children_edges

  def add_insplit_range_edges(self, r_inv_e, new_vertices, edges, v):
    '''
    Function for edges whose range is v   ( r(e) = v  or e in r^{-1}(v) )
    Must partition edges into nonempty sets E1, E2
    '''
    # edge in Ei need their range set to v^i
    for e in edges:
      if e in self.E1:
        e.r = new_vertices[0]
      if e in self.E2:
        e.r = new_vertices[1]

    return edges
