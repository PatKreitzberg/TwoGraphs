import sys

from TwoGraph import *
from calculate_h1 import calculate_h1_gemini, calculate_h1_claude

def add_insplit_edges(g_edges, s_inv_e, v, new_vertices):
  edges = set(g_edges) # list
  edges -= s_inv_e # set
  child_edge_function= {e:[e] for e in edges} # these edges are not modified at all
  new_edges = set()

  print("WARNING: Do we redirect edges which now point to v^1, v^2 ????")

  # ADD NEW SOURCE EDGES
  for e in s_inv_e:
    child_edge_function[e] = list()
    for j in [0,1]:
      label = e.label+"^" + str(j)
      s = new_vertices[j]
      r = e.r
      new_e = Edge(label, s, r)
      new_e.degree_index = e.degree_index

      edges.add(new_e)
      new_edges.add(new_e)
      child_edge_function[e].append(new_e)

  return edges, child_edge_function, new_edges

def add_insplit_commuting_squares(v, old_commuting_squares, child_edge_function, new_edges):
  # we have af ~ eb if
  # 1. The parent squares in the original graph commute
  # 2. The range(af) = range(eb)
  # 3. The source(af) = source(eb)


  # Commuting squares which do not have edges whose source is v in
  # original graph. Insplitting preserves range of commuting squares
  # and if the commuting square does not contain an edge with source
  # of v then the commuting square is unphased. Just the ede vertex is
  # changed.
  new_commuting_squares = set()
  for cs in old_commuting_squares:
    include = True
    for new_edge in new_edges:
      if new_edge in cs:
        include = False
        continue
    if include:
      new_commuting_squares.add(cs)

  commuting_squares_to_inspect = set().union(*[new_edge.commuting_squares for new_edge in new_edges])

  for cs in commuting_squares_to_inspect:
    (a,e),(f,b) = cs.lhs, cs.rhs
    for a_ in child_edge_function[a]:
        for e_ in child_edge_function[e]:
            for f_ in child_edge_function[f]:
                for b_ in child_edge_function[b]:
                  if check_commuting_square_is_valid(a_, e_, f_, b_):
                    new_commuting_squares.add(CommutingSquare(a_, e_, f_, b_))

  return list(new_commuting_squares)

def check_commuting_square_is_valid(a, e, f, b):
  # a e ~ f b
  # need:
  # 1. We have two paths
  #     range of e is source of a
  #     range of b is source of f
  # 2. The two paths have same source and range:
  #     source of e is source of b
  #     range of a is range of f
  return (a.s == e.r) and (a.r == f.r) and (e.s == b.s) and (b.r == f.s)


def add_insplit_vertices(v, vertices, new_vertices):
  vertices = set(vertices)
  vertices.remove(v)
  vertices |= set(new_vertices)
  return vertices

def insplit(g, v):
  # g is a TwoGraph
  # v is a str

  # When we create this new graph order of edges and vertices don't matter

  # ADD NEW VERTICES
  new_vertices = [v+'^1', v+'^2']
  vertices_as_set = add_insplit_vertices(v, g.vertices, new_vertices)

  # MODIFY EDGES

  # ADD NEW RANGE EDGES (r(e) = v)
  # Partitions r_inv_e
  r_inv_e = set(g.range_inverse_of_vertex(v))
  E1, E2 = partition_pairs(r_inv_e)
  assert len(E1) > 0
  assert len(E2) > 0
  for e in E1:
    e.r = new_vertices[0]
  for e in E2:
    e.r = new_vertices[1]


  # add new source edges (s(e) = v)
  s_inv_e = set(g.source_inverse_of_vertex(v))
  edges, child_edge_function, new_edges = add_insplit_edges(g.edges, s_inv_e, v, new_vertices)


  # COMMUTING SQUARES
  commuting_squares = add_insplit_commuting_squares(v, g.commuting_squares, child_edge_function, new_edges)

  g_insplit = TwoGraph({
    'vertices':list(vertices_as_set),
    'edge_label_to_edge':{edge.label:edge for edge in edges},
    'commuting_squares':commuting_squares
  })
  return g_insplit



def partition_pairs(r_inv_e):
  commuting_squares_with_range_v = set()
  for e in r_inv_e:
    commuting_squares_with_range_v |= e.range_of_commuting_squares

  S = set()
  for cs in commuting_squares_with_range_v:
    S.add( tuple(sorted(cs.r)) )


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


def calc_homology_and_insplit_homology(path, plot):
  g = TwoGraph(path)

  print("######## INSPLITTING ###########")

  g_i = insplit(g, 'v')
  g_i.draw_graph()
  return exit()



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

  if plot:
    g.draw_graph()

if __name__ == "__main__":
  print()
  if len(sys.argv) < 2:
    print("Usage: python two_graph.py <file> <file2>....")
    sys.exit(1)

  plot = False
  if 'plot' in sys.argv:
   plot=True


  for path in sys.argv[1:]:
    if path != 'plot':
      calc_homology_and_insplit_homology(path, plot)
