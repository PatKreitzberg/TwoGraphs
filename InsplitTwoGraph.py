from TwoGraph import TwoGraph
from Edge import Edge
from CommutingSquare import CommutingSquare

class InsplitTwoGraph(TwoGraph):
  def __init__(self, g, v, E1=None, E2=None):
    # og_graph is the graph we are going to insplit
    # v is the vertex at which we insplit

    #Add vertices to graph
    new_vertices = [str(v)+'^1', str(v)+'^2']
    vertices_as_set = self.add_insplit_vertices(v, g.vertices, new_vertices)
    self.v_as_s = vertices_as_set

    # MODIFY EDGES
    s_inv_e = set(g.source_inverse_of_vertex(v))
    edges, edge_to_children_edges = self.add_insplit_edges(set(g.range_inverse_of_vertex(v)), s_inv_e, set(g.edges), new_vertices, v, E1, E2)

    # COMMUTING SQUARES
    commuting_squares = self.add_insplit_commuting_squares(v, g.commuting_squares, edge_to_children_edges)

    load_from = {
      'vertices':list(vertices_as_set),
      'edge_label_to_edge':{edge.label:edge for edge in edges},
      'commuting_squares':commuting_squares
    }
    super().__init__(load_from)


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

  def add_insplit_edges(self, r_inv_e, s_inv_e, edges, new_vertices, v, E1, E2):
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

  def add_insplit_range_edges(self, r_inv_e, new_vertices, edges, v, E1, E2):
    '''
    Function for edges whose range is v   ( r(e) = v  or e in r^{-1}(v) )
    Must partition edges into nonempty sets E1, E2
    '''
    if (E1 is None) or (E2 is None):
      E1, E2 = self.partition_pairs(r_inv_e)

    assert len(E1&E2) == 0
    assert len(E1) > 0
    assert len(E2) > 0

    # edge in Ei need their range set to v^i
    for e in edges:
      if e in E1:
        e.r = new_vertices[0]
      if e in E2:
        e.r = new_vertices[1]

    return edges


  def partition_pairs(self, r_inv_e):
    commuting_squares_with_range_v = set()
    for e in r_inv_e:
      commuting_squares_with_range_v |= e.range_of_commuting_squares

    S = set()
    for cs in commuting_squares_with_range_v:
      S.add( tuple(sorted(cs.r)) )
    print("S", S)

    # 1. Map every element to its "parent" (initially itself)
    parent = {}

    def find(i):
      if parent[i] == i:
        return i
      parent[i] = find(parent[i]) # Path compression
      return parent[i]

    def union(i, j):
      root_i = find(i)
      root_j = find(j)
      if root_i != root_j:
        parent[root_i] = root_j

    # Get all unique elements
    elements = set().union(*S)
    for el in elements:
      parent[el] = el

    # 2. Union the elements in each pair
    for pair in S:
      p_list = list(pair)
      # Since it's a pair, we union the first element with the second
      union(p_list[0], p_list[1])

    # 3. Group elements by their root parent
    components = {}
    for el in elements:
      root = find(el)
      if root not in components:
        components[root] = set()
      components[root].add(el)

    # 4. Partition into E1 and E2
    # We take the first component found for E1,
    # and put everything else into E2.
    groups = list(components.values())

    if not groups:
      return set(), set()

    E1 = groups[0]
    E2 = set().union(*groups[1:]) if len(groups) > 1 else set()

    # Put any elements not in the partitions in a partition
    r_inv_e -= E1
    r_inv_e -= E2
    if len(r_inv_e) > 0:
      i = 0
      for e in r_inv_e:
        if (i%2)== 0:
          E1.add(e)
        else:
          E2.add(e)
        i += 1

    return E1, E2
