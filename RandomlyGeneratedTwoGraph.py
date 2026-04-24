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
    R_degree = 1
    B_degree = 2
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

    commuting_squares = self.get_commuting_squares(R_path_matrix, B_path_matrix)

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


  def partition_for_insplit(self, R,B,range_vertex_to_commuting_square):
    '''
    v is a vertex
    Need to find all edges with range v
    Partition them then create commuting squares

    Need at least four source edges! Then can always insplit. So need a vertex that has degree >= 4
    '''
    print("R")
    print(R.adj_matrix)
    print()
    print("B")
    print(B.adj_matrix)
    print()

   # Finding a vertex with at least four in-edges
    for v in range(self.n):
      total_degree = 0
      for row in range(self.n):
        total_degree += len(R[row][v])
        total_degree += len(B[row][v])
      if total_degree >= 4:
        self.can_insplit = True
        self.insplit_v = v
        print(f"At vertex {self.insplit_v} degree {total_degree}")
        break

    # If we can't insplit don't bother
    if not self.can_insplit:
      return None

    degree_1_source_edge_to_paths = dd(set)
    degree_2_source_edge_to_paths = dd(set)
    # Need to partition paths by source edge and the degree of that source edge

    for cs in range_vertex_to_commuting_square[self.insplit_v]:
      s1,r1 = cs.path1
      s2,r2 = cs.path2

      if s1.degree == 1:
        degree_1_source_edge_to_paths[s1].add(cs.path1)
        degree_2_source_edge_to_paths[s2].add(cs.path2)
      else:
        degree_2_source_edge_to_paths[s1].add(cs.path1)
        degree_1_source_edge_to_paths[s2].add(cs.path2)

    X = [(len(paths), se) for se,paths in degree_1_source_edge_to_paths.items()]
    Y = [(len(paths), se) for se,paths in degree_2_source_edge_to_paths.items()]

    equal_subset = find_first_equal_subset(X,Y)
    sub_X = equal_subset['subset_A']
    sub_Y = equal_subset['subset_B']

    all_commuting_squares = set()
    # sub_X and sub_Y will create a partition of source edges that allow us to insplit
    all_paths_degree_1 = set.union(*[degree_1_source_edge_to_paths[se] for _,se in sub_X])
    all_paths_degree_2 = set.union(*[degree_2_source_edge_to_paths[se] for _,se in sub_Y])
    # these commuting squares guarantees we can insplit
    for p1,p2 in zip(all_paths_degree_1, all_paths_degree_2):
      all_commuting_squares.add(CommutingSquare(*p1, *p2))

    # Add all other commuting squares
    for range_vertex, commuting_squares in range_vertex_to_commuting_square.items():
      if range_vertex == self.insplit_v:
        continue
      all_commuting_squares |= commuting_squares

    return all_commuting_squares

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
          print('CHECK Source and Range are correct order RB path', CommutingSquare(r1,s1, r2, s2))

          range_vertex = r1.r
          range_vertex_to_commuting_square[range_vertex].add(CommutingSquare(r1,s1, r2, s2))

    return self.partition_for_insplit(R, B, range_vertex_to_commuting_square)
