import numpy as np
from collections import defaultdict as dd
import multiprocessing
from functools import partial

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

    #R,B = self.gen()
    R,B = None,None

    # no torsion but all homology differs
    if False:
      self.n = 3
      self.z = 3
      R = [[1, 0, 2],[0, 1, 0],[0, 2, 1]]
      B = [[2, 4, 4], [0, 2, 0],[0, 4, 2]]
    elif False:
      # torsion! and H1, H2 differ
      self.n = 5
      self.z = 3
      R = [[1, 2, 0, 0, 0],[0, 1, 0, 0, 0],[0, 0, 1, 0, 0],[0, 0, 0, 1, 0],[2, 1, 0, 0, 1]]
      B = [[2, 4, 0, 0, 0],[0, 2, 0, 0, 0],[0, 0, 2, 0, 0],[0, 0, 0, 2, 0],[4, 6, 0, 0, 2]]
    elif 0:
      self.n = 2
      self.z = 3
      R = [[1,1],[1, 1]]
      B = R
    elif 0:
      B = R
    elif 0:
      R = [[2,0],[0, 0]]
      B = R
    elif 1:
      R = [[2]]
      B = R





    assert R is not None
    assert B is not None

    # Setup matrices
    self.R_path_matrix = PathMatrix(R, self.R_degree)
    self.B_path_matrix = PathMatrix(B, self.B_degree)
    self.RB_paths_matrix = self.R_path_matrix*self.B_path_matrix
    self.BR_paths_matrix = self.B_path_matrix*self.R_path_matrix

    self.print_path_matrix(self.R_path_matrix.adj_matrix, "Red edge adjacency matrix")
    self.print_path_matrix(self.B_path_matrix.adj_matrix, "Blue edge adjacency matrix")

    if type(self) is RandomlyGeneratedTwoGraph:
      self.commuting_squares = self.get_commuting_squares()
      assert self.commuting_squares is not None

      n_commuting_squares = len(self.commuting_squares)
      n_rb_paths = sum([sum([len(col) for col in row]) for row in self.BR_paths_matrix])
      assert n_commuting_squares == n_rb_paths

      self.vertices = [i for i in range(n)]
      self.edges = {edge for edge in self.R_path_matrix.edges + self.B_path_matrix.edges}
      self.calculate_boundary_matrices()


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

  def commuting_square_for_vertices(self, RB_path, BR_path, commuting_squares):
    s1,r1 = RB_path
    s2,r2 = BR_path

    assert s1.s == s2.s
    assert r1.r == r2.r
    commuting_squares.append(CommutingSquare(r1,s1, r2, s2))
    return commuting_squares

  def commuting_squares_for_vertices(self, source_v, range_v, commuting_squares):
    assert len(self.RB_paths_matrix[source_v][range_v]) == len(self.BR_paths_matrix[source_v][range_v])
    for RB_path, BR_path  in zip(self.RB_paths_matrix[source_v][range_v], self.BR_paths_matrix[source_v][range_v]):
      commuting_squares = self.commuting_square_for_vertices(RB_path, BR_path, commuting_squares)
    return commuting_squares

  def get_commuting_squares(self):
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
    for source_v in range(self.n): # source vertex
      for range_v in range(self.n): # range vertex
        commuting_squares = self.commuting_squares_for_vertices(source_v, range_v, commuting_squares)
    return commuting_squares


  def gen(self, attempt=0):
    """
      Generates a unimodular matrix A and a commuting matrix B.
      """
    # 1. Create a starting triangular matrix with 1 or -1 on diagonal
    # This guarantees det(A) = 1 or -1
    A = np.eye(self.n, dtype=object)
    diagonals = np.random.choice([1, -1], size=self.n)
    for i in range(self.n):
      A[i, i] = diagonals[i]

      # 2. Apply random shear operations (adding one row to another)
      # These operations preserve the determinant.
      for _ in range(self.n * 2):
        i, j = np.random.choice(self.n, size=2, replace=False)
        # factor is what actually changes the values
        factor = np.random.randint(1, 3) # Small factors to keep entries manageable
        A[i] += factor * A[j]

      # 3. Ensure entries are within roughly [0, z]
      # (Strictly bound constraints are hard with det=1, but we can take modulo or adjust)
      A = A.astype(int) % (self.z + 1)

      # Re-verify/Force determinant for the example (simple construction)
      # Let's use a simpler approach for a guaranteed clean result:
      # A = I + Nilpotent or just a small shear sequence.

      # --- Practical Approach for AB = BA ---
      # B = cI is trivial. Let's make B = A^2 + kI to ensure positivity.
      # We add a large enough constant k to make all entries positive.
      B_raw = np.dot(A, A)
      min_val = np.min(B_raw)
      k = abs(min_val) + 1 if min_val <= 0 else 1
      B = B_raw + k * np.eye(self.n, dtype=int)
      if (np.array_equal(np.dot(A, B), np.dot(B, A))):
        return A, B
      else:
        return self.gen(attempt+1)

  def gen_circulant_matrices(self):
    # THESE SEEM TO NEVER BE ABLE TO INSPLIT
    try:
      # Get user input
      n = self.n
      z = self.z

      if n <= 0 or z < 0:
        print("Please enter a positive size and a non-negative bound.")
        return

      # A circulant matrix is defined by its first row.
      # Any two circulant matrices A and B will satisfy AB = BA.

      # Generate random first rows with entries in [0, z]
      row_a = np.random.randint(0, z + 1, size=n)
      row_b = np.random.randint(0, z + 1, size=n)

      # Construct matrices A and B by shifting the first row
      A = np.array([np.roll(row_a, i) for i in range(n)])
      B = np.array([np.roll(row_b, i) for i in range(n)])

      print("\nMatrix A:")
      print(A)

      print("\nMatrix B:")
      print(B)

      # Verification
      AB = np.dot(A, B)
      BA = np.dot(B, A)

      print("\nVerification (AB == BA):")
      if np.array_equal(AB, BA):
        print("Success! AB is equal to BA.")
      else:
        print("Failure. The matrices do not commute.")

      # Optional: Print the product result
      # print("\nProduct AB:")
      # print(AB)
      return A,B
    except ValueError:
      print("Invalid input. Please enter integers only.")
