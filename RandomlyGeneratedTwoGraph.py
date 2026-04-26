import numpy as np
import itertools
import collections
from collections import defaultdict as dd

from find_first_equal_subset import find_first_equal_subset
from TwoGraph import TwoGraph
from PathMatrix import PathMatrix
from CommutingSquare import CommutingSquare

class RandomlyGeneratedTwoGraph(TwoGraph):
  def __init__(self, n, z):
    super().__init__()
    # n is number of vertices
    # z is an upper limit on number of edges in graph between any two (maybe non-distinct) vertices
    self.n = n
    self.z = z
    self.is_legit = True
    attempts = 0

    R,B = self.generate_adjacency_matrices()
    assert R is not None
    assert B is not None

    # Setup matrices
    self.R_path_matrix = PathMatrix(R, self.R_degree)
    self.B_path_matrix = PathMatrix(B, self.B_degree)
    self.RB_paths_matrix = self.R_path_matrix*self.B_path_matrix
    self.BR_paths_matrix = self.B_path_matrix*self.R_path_matrix

    # Find vertex which can be partitioned at
    v,E1,E2 = self.find_insplit_vertex()
    commuting_squares = self.get_commuting_squares(v, E1,E2)
    assert commuting_squares is not None
    assert v is not None
    assert len(E1) > 0
    assert len(E2) > 0

    n_commuting_squares = len(commuting_squares)
    n_rb_paths = sum([sum([len(col) for col in row]) for row in self.BR_paths_matrix])
    assert n_commuting_squares == n_rb_paths

    self.vertices = {i for i in range(n)}
    self.edges = {edge for edge in self.R_path_matrix.edges + self.B_path_matrix.edges}
    self.commuting_squares = commuting_squares
    self.calculate_boundary_matrices()

  def find_insplit_vertex(self):
    for v in range(self.n):
      E1, E2 =  self.matching_partition(v)
      if len(E1)>0 and len(E2)>0:
        return v, E1, E2
    return None,{},{}

  def generate_adjacency_matrices(self, attempt=0):
    ''''
    Generates adjacency matrices R and B
    R is just for red edges, B for blue edges

    CONSTRAINTS:
    RB = BR so that for every length two red-blue path there is an equivalent blue-red path
    Elements of R, B are nonnegative integers between 0 and z, roughly. For the first generated matrix this is
    always true; for the scond matrix it is harder to guarantee because it must be tailored to commute with
    the first matrix
    '''
    if attempt > 5:
      return None,None

    # 1. Generate A with det = 1 or -1
    # We use a lower triangular matrix with 1s on diagonal and small off-diagonals
    A = np.eye(self.n, dtype=int)

    # Fill lower triangle with small integers to keep growth controlled
    np.random.seed(0)
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
    if sum(sum(np.dot(A, B) - np.dot(A, B))) != 0:
      return self.generate_adjacency_matrices(attempt=attempt+1)
    return A, B

  def get_commuting_squares(self, v, E1, E2):
    '''
    Inputs:
    R,B PathMatrix objects
    They are n x n matrices
      R[i][j] is a list of paths of length 1 from v_i to v_j
      path of the form
    '''
    # Of the form ('r', i, j), ('b', u, v)  with
    # source edge being ('b',u,v)
    # range edge being ('r',i,j)
    commuting_squares = []
    range_vertex_to_commuting_square = dd(set)
    for row in range(self.n): # source vertex
      for col in range(self.n): # range vertex
        if col == v:
          continue

        assert len(self.RB_paths_matrix[row][col]) == len(self.BR_paths_matrix[row][col])

        for RB_path, BR_path  in zip(self.RB_paths_matrix[row][col], self.BR_paths_matrix[row][col]):
          s1,r1 = RB_path
          s2,r2 = BR_path

          assert s1.s == s2.s
          assert r1.r == r2.r
          commuting_squares.append(CommutingSquare(r1,s1, r2, s2))

    # Commuting squares for the insplit vertex
    E1_RB_paths,E1_BR_paths = set(),set()
    E2_RB_paths,E2_BR_paths = set(),set()

    for row in range(self.n):
      E1_RB_paths, E2_RB_paths = self.partition_paths(E1, E2, self.RB_paths_matrix[row][v])
      E1_BR_paths, E2_BR_paths = self.partition_paths(E1, E2, self.BR_paths_matrix[row][v])

      for E_RB_paths,E_BR_paths in [(E1_RB_paths, E1_BR_paths), (E2_RB_paths, E2_BR_paths)]:
        for RB_path, BR_path in zip(E_RB_paths,E_BR_paths):
            s1,r1 = RB_path
            s2,r2 = BR_path
            assert s1.s == s2.s
            assert r1.r == r2.r
            commuting_squares.append(CommutingSquare(r1,s1, r2, s2))

    return commuting_squares

  def partition_paths(self, E1, E2, paths):
    E1_paths, E2_paths = set(), set()
    for se,re in paths:
      if re in E1:
        E1_paths.add((se,re))
      elif  re in E2:
        E2_paths.add((se,re))
      else:
        assert False
    return E1_paths, E2_paths


  def get_source_profile_vectors(self, v):
    '''
    Want to get vectors Ve. One for each edge e in r^{-1}(v). They will have size 1xn.
    Ve[i] = number of paths from vertex i to the source vertex of e such that the path is a red-blue path
    (assuming e is blue edge)
    '''

    # Edges whose range is v
    incoming_red_edges  = [e for e in self.R_path_matrix.edges if e.r == v]
    for re in incoming_red_edges:
      assert re.degree == self.R_degree
    red_range_edge_to_num_paths_to_ui  = {e:[0]*self.n for e in incoming_red_edges}

    incoming_blue_edges = [e for e in self.B_path_matrix.edges if e.r == v]
    for be in incoming_blue_edges:
      assert be.degree == self.B_degree
    blue_range_edge_to_num_paths_to_ui = {e:[0]*self.n for e in incoming_blue_edges}

    for i in range(self.n):
      # For each vertex i find the number of blue-red paths which go
      # from vertex i to v through blue edge f
      for _,blue_range_edge in self.RB_paths_matrix[i][v]:
        assert blue_range_edge.degree == self.B_degree
        blue_range_edge_to_num_paths_to_ui[blue_range_edge][i] += 1

      for _,red_range_edge in self.BR_paths_matrix[i][v]:
        assert red_range_edge.degree == self.R_degree
        red_range_edge_to_num_paths_to_ui[red_range_edge][i] += 1

    blue_range_edge_to_num_paths_to_ui = {e: tuple(vec) for e,vec in blue_range_edge_to_num_paths_to_ui.items()}
    red_range_edge_to_num_paths_to_ui = {e: tuple(vec) for e,vec in red_range_edge_to_num_paths_to_ui.items()}
    return incoming_red_edges, incoming_blue_edges, red_range_edge_to_num_paths_to_ui, blue_range_edge_to_num_paths_to_ui


  def matching_partition(self, v):
    '''Get source profile vectors then do a subset sum problem where
    we want a subset of the red edge to number of paths vectors and a
    subset of the blue edge to number of path vectors such all the
    blue sum to the same vector as all the red. If two proper subsets
    exist then they are contain the edges that should be in E1. E2 is
    the rest of the edges.
    '''

    R_edges, B_edges, R_vecs, B_vecs = self.get_source_profile_vectors(v)
    all_incoming = list(set(R_edges) | set(B_edges))

    if len(all_incoming) < 2:
      return None

    # Shortcut: Check for edges that participate in zero commuting squares
    for edge in all_incoming:
      vec = R_vecs.get(edge) or B_vecs.get(edge)
      if all(x == 0 for x in vec):
        E1 = {edge}
        E2 = set(all_incoming) - E1
        return E1, E2

    # Vector Subset Sum Search
    m, n = len(R_edges), len(B_edges)

    # We look for a subset of red edges I and blue edges J
    # such that the sum of their requirement vectors is identical.
    for r_size in range(m + 1):
      for b_size in range(n + 1):
        if (r_size == 0 and b_size == 0) or (r_size == m and b_size == n):
          continue

        # Optimization: Pre-calculate combinations of Red subsets
        for I_sub in itertools.combinations(range(m), r_size):
          red_sum = [0] * self.n
          for idx in I_sub:
            for k in range(self.n):
              red_sum[k] += R_vecs[R_edges[idx]][k]

          # Check against all combinations of Blue subsets
          for J_sub in itertools.combinations(range(n), b_size):
            blue_sum = [0] * self.n
            for idx in J_sub:
              for k in range(self.n):
                blue_sum[k] += B_vecs[B_edges[idx]][k]

            if red_sum == blue_sum:
              E1 = set(R_edges[i] for i in I_sub) | set(B_edges[j] for j in J_sub)
              E2 = set(all_incoming) - E1
              return E1, E2

    return {},{}
