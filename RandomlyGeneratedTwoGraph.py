import numpy as np
from collections import defaultdict as dd

from find_first_equal_subset import find_first_equal_subset
from TwoGraph import TwoGraph
from PathMatrix import PathMatrix
from CommutingSquare import CommutingSquare

class RandomlyGeneratedTwoGraph(TwoGraph):
  def __init__(self, n, z):
    # n is number of vertices
    # z is an upper limit on number of edges in graph between any two (maybe non-distinct) vertices
    self.can_insplit = False
    self.n = n
    self.z = z
    self.R_degree = 1
    self.B_degree = 2
    attempts = 0

    while True:
      # Make sure R, B truly commute.
      R,B = self.generate_adjacency_matrices()
      if sum(sum(np.dot(R, B) - np.dot(B, R))) == 0:
        break
      attempts += 1
      if attempts > 5:
        assert False

    R_path_matrix = PathMatrix(R, R_degree)
    B_path_matrix = PathMatrix(B, B_degree)

    commuting_squares, E1, E2 = self.get_commuting_squares(R_path_matrix, B_path_matrix)
    self.E1 = E1
    self.E2 = E2

    if commuting_squares is None:
      return

    print("Okay creating TwoGraph from Random graph now")
    load_from = {
      'vertices':[i for i in range(n)],   # Vertices are just basic integers 0,1,...,n-1
      'edge_label_to_edge':{edge.label:edge for edge in R_path_matrix.edges + B_path_matrix.edges},
      'commuting_squares':commuting_squares
    }
    super().__init__(load_from)




  def generate_adjacency_matrices(self):
    ''''
    Generates adjacency matrices R and B
    R is just for red edges, B for blue edges

    CONSTRAINTS:
    RB = BR so that for every length two red-blue path there is an equivalent blue-red path
    Elements of R, B are nonnegative integers between 0 and z, roughly. For the first generated matrix this is
    always true; for the scond matrix it is harder to guarantee because it must be tailored to commute with
    the first matrix
    '''
    # 1. Generate A with det = 1 or -1
    # We use a lower triangular matrix with 1s on diagonal and small off-diagonals
    A = np.eye(self.n, dtype=int)

    # Fill lower triangle with small integers to keep growth controlled
    for i in range(self.n):
        for j in range(i):
            A[i, j] = np.random.randint(0, self.z // 2)

    # Apply a few row/column swaps or additions to "hide" the triangular structure
    # This keeps det(A) as 1 * 1 * ... * 1 = 1
    for _ in range(self.n):
        idx = np.random.permutation(self.n)
        A[idx[0]] += A[idx[1]]

    # 2. Generate B = A + kI
    # Find the minimum value in A to determine the smallest possible k
    min_a = np.min(A)
    # k must make all entries > 0. If min_a is 0, k=1. If min_a is -5, k=6.
    k = max(1, 1 - min_a)

    B = A + k * np.eye(self.n, dtype=int)
    return A, B


  def find_vertex_with_incoming_degree_at_least_four(self, R, B):
    # Finding a vertex with at least four in-edges
    for range_vertex in range(self.n):
      total_edges_into_range_vertex = 0
      for source_vertex in range(self.n):
        total_edges_into_range_vertex += len(R[source_vertex][range_vertex]) + len(B[source_vertex][range_vertex])
      if total_edges_into_range_vertex >= 4:
        self.can_insplit = True
        return range_vertex
    return None

  def partition_for_insplit(self, R,B,range_vertex_to_commuting_square):
    '''
    v is a vertex
    Need to find all edges with range v
    Partition them then create commuting squares

    Need at least four source edges! Then can always insplit. So need a vertex that has degree >= 4

    Procedure:
    Take all edges with range r
    Separate into degree 1 and degree 2
    Need to get partitions E1, E2 such that
    1. E1 and E2 have as many degree 1 edges as degree 2 edges
    2. E1, E2 are nonempty, no shared elements, partittion range edges of v

    Make commuting edges where r1 s1 ~ r2 s2 => r1, r2 both in E1 or both in E2 !
    '''
    self.insplit_v = self.find_vertex_with_incoming_degree_at_least_four(R,B)

    # If we can't insplit don't bother
    if not self.can_insplit:
      return None

    degree_1_range_edge_to_paths = dd(set)
    degree_2_range_edge_to_paths = dd(set)
    # Need to partition paths by range edge and the degree of that range edge

    for cs in range_vertex_to_commuting_square[self.insplit_v]:
      r1, s1 = cs.path1
      r2, s2 = cs.path2

      if r1.degree == 1:
        degree_1_range_edge_to_paths[r1].add(cs.path1)
        degree_2_range_edge_to_paths[r2].add(cs.path2)
      else:
        degree_2_range_edge_to_paths[r1].add(cs.path1)
        degree_1_range_edge_to_paths[r2].add(cs.path2)

    X = [(len(paths), re) for re,paths in degree_1_range_edge_to_paths.items()]
    Y = [(len(paths), re) for re,paths in degree_2_range_edge_to_paths.items()]

    equal_subset = find_first_equal_subset(X,Y)
    E1_degree_1_edges = [edge for _,edge in equal_subset['subset_A']]
    E1_degree_2_edges = [edge for _,edge in equal_subset['subset_B']]

    E2_degree_1_edges = [edge for _,edge in X if edge not in E1_degree_1_edges]
    E2_degree_2_edges = [edge for _,edge in Y if edge not in E1_degree_2_edges]
    return (E1_degree_1_edges, E1_degree_2_edges), (E2_degree_1_edges, E2_degree_2_edges)

  def get_commuting_squares(self, R, B):
    '''
    Inputs:
    R,B PathMatrix objects
    They are n x n matrices
      R[i][j] is a list of paths of length 1 from v_i to v_j
      path of the form
    '''
    print("WARNING: enforce we can insplit at vertex")
    print("WARNING: Change to sets instead of lists otherwise commuting squares will always be boring I think?")

    # Of the form ('r', i, j), ('b', u, v)  with
    # source edge being ('b',u,v)
    # range edge being ('r',i,j)
    RB_paths_matrix = R*B
    BR_paths_matrix = B*R
    commuting_squares = []

    range_vertex_to_commuting_square = dd(set)

    for row in range(self.n): # source vertex
      for col in range(self.n): # range vertex
        assert len(RB_paths_matrix[row][col]) == len(BR_paths_matrix[row][col])
        for RB_path, BR_path  in zip(RB_paths_matrix[row][col], BR_paths_matrix[row][col]):
          s1,r1 = RB_path
          s2,r2 = BR_path

          assert s1.s == s2.s
          assert r1.r == r2.r

          range_vertex_to_commuting_square[r1.r].add(CommutingSquare(r1,s1, r2, s2))

    (E1_red_edges, E1_blue_edges), (E2_red_edges, E2_blue_edges) = self.partition_for_insplit(R, B, range_vertex_to_commuting_square)

    all_commuting_squares = set()
    all_commuting_squares |= self.get_commuting_squares_for_partition(E1_red_edges, E1_blue_edges, R, B)
    all_commuting_squares |= self.get_commuting_squares_for_partition(E2_red_edges, E2_blue_edges, R, B)

    for range_vertex, commuting_squares in range_vertex_to_commuting_square.items():
      if range_vertex == self.insplit_v:
        continue
      all_commuting_squares |= commuting_squares

    E1 = E1_red_edges + E1_blue_edges
    E2 = E2_red_edges + E2_blue_edges
    return commuting_squares, E1, E2


  def get_commuting_squares_for_partition(self, E_red_edges, E_blue_edges, R, B):
    red_blue_paths = set()
    for red_edge in E_red_edges: # Red edges
      # it is the range edge so we need to look at the other matrix to get source edge
      for source_vertex in len(B.n):
        for blue_edge in B[source_vertex][red_edge.s]:
          red_blue_paths.add((red_edge, blue_edge))

    blue_red_paths = set()
    for blue_edge in E_blue_edges:
      for source_vertex in len(R.n):
        for red_edge in R[source_vertex][blue_edge.s]:
          blue_red_paths.add((blue_edge, red_edge))
    commuting_squares = set()

    assert len(red_blue_paths) == len(blue_red_paths)
    for (r1,s1), (r2,s2) in zip(red_blue_paths, blue_red_paths):
      commuting_squares.add(CommutingSquare(r1,s1,r2,s2))
    return commuting_squares
