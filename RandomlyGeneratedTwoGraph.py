import numpy as np
import itertools
import collections
from collections import defaultdict as dd
import multiprocessing
from functools import partial

from find_first_equal_subset import find_first_equal_subset
from TwoGraph import TwoGraph
from PathMatrix import PathMatrix
from CommutingSquare import CommutingSquare

def _worker_wrapper(v, obj):
  """
  Top-level helper function to call the method on the object.
  This avoids pickling issues with bound methods.
  """
  E1, E2 = obj.matching_partition_claude(v)
  print("V=",v,"done")
  if len(E1) > 0 and len(E2) > 0:
    return v, E1, E2
  return None

class RandomlyGeneratedTwoGraph(TwoGraph):
  def __init__(self, n, z):
    super().__init__()
    # n is number of vertices
    # z is an upper limit on number of edges in graph between any two (maybe non-distinct) vertices
    self.n = n
    self.z = z
    self.is_legit = True

    #R,B = self.generate_adjacency_matrices()
    R,B = self.gen()

    # WITH TORSION
   # R,B with H1, H2 different. Genereated with the gen() function, not generate_adjacency_matrices() function
    # R = [[1, 2, 0, 0, 0],[0, 1, 0, 0, 0],[0, 0, 1, 0, 0],[0, 0, 0, 1, 0],[2, 1, 0, 0, 1]]
    # B = [[2, 4, 0, 0, 0],[0, 2, 0, 0, 0],[0, 0, 2, 0, 0],[0, 0, 0, 2, 0],[4, 6, 0, 0, 2]]

    # WITHOUT TORSION
    # R,B with H1, H2 different. Genereated with the gen() function, not generate_adjacency_matrices() function
    # R = [[1, 0, 0, 1],[2, 1, 0, 0],[0, 0, 1, 1],[0, 0, 0, 1]]
    # B = [[2, 0, 0, 2],[4, 2, 0, 2],[0, 0, 2, 2],[0, 0, 0, 2]]

    assert R is not None
    assert B is not None

    # Setup matrices
    self.R_path_matrix = PathMatrix(R, self.R_degree)
    self.B_path_matrix = PathMatrix(B, self.B_degree)
    self.RB_paths_matrix = self.R_path_matrix*self.B_path_matrix
    self.BR_paths_matrix = self.B_path_matrix*self.R_path_matrix

    self.print_path_matrix(self.R_path_matrix.adj_matrix, "Red edge adjacency matrix")
    self.print_path_matrix(self.B_path_matrix.adj_matrix, "Blue edge adjacency matrix")

    # Find vertex which can be partitioned at
    #v,E1,E2 = self.find_insplit_vertex()
    v,E1,E2 = self.find_insplit_vertex_opt()
    print(f"v is {v}")
    assert v is not None
    assert len(E1) > 0
    assert len(E2) > 0

    commuting_squares = self.get_commuting_squares(v, E1,E2)
    assert commuting_squares is not None

    n_commuting_squares = len(commuting_squares)
    n_rb_paths = sum([sum([len(col) for col in row]) for row in self.BR_paths_matrix])
    assert n_commuting_squares == n_rb_paths

    self.v = v
    self.E1 = E1
    self.E2 = E2
    self.vertices = {i for i in range(n)}
    self.edges = {edge for edge in self.R_path_matrix.edges + self.B_path_matrix.edges}
    self.commuting_squares = commuting_squares
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

    # Edges whose range is v
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
      print("Taking the long way...", r_size)
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

  def find_insplit_vertex(self):
    # Determine the number of available CPU cores
    num_cores = multiprocessing.cpu_count()

    # Use a partial to pass 'self' into the worker
    worker = partial(_worker_wrapper, obj=self)

    with multiprocessing.Pool(processes=num_cores) as pool:
      # imap_unordered is faster and allows us to stop early
      for result in pool.imap_unordered(worker, range(self.n)):
        if result is not None:
          # Terminate the pool immediately once we find a match
          pool.terminate()
          return result
    return None, {}, {}

  def find_insplit_vertex_opt(self):
    v_to_graph_degree = {v:0 for v in range(self.n)}
    for se in range(self.n):
      for re in range(self.n):
        v_to_graph_degree[re] += self.R_path_matrix.adj_matrix[se][re]
        v_to_graph_degree[re] += self.B_path_matrix.adj_matrix[se][re]

    v_sorted_by_degree = sorted( [(int(d),v) for v,d in v_to_graph_degree.items()] )

    for _,v in v_sorted_by_degree:
      print("Trying with v=",v)
      E1, E2 =  self.matching_partition_opt(v)
      if len(E1)>0 and len(E2)>0:
        return v, E1, E2
    return None,{},{}



  def matching_partition_opt(self, v):
      R_edges, B_edges, R_vecs, B_vecs = self.get_source_profile_vectors(v)
      all_incoming = list(set(R_edges) | set(B_edges))

      if len(all_incoming) < 2:
          return None

      # 1. Zero-vector shortcut (remains high priority)
      for edge in all_incoming:
          vec = R_vecs.get(edge) or B_vecs.get(edge)
          if all(x == 0 for x in vec):
              E1 = {edge}
              return E1, set(all_incoming) - E1

      m, n = len(R_edges), len(B_edges)

      # 2. Build a Map of all possible Red subset sums
      # Key: Sum tuple, Value: Set of edge indices
      red_sum_map = {}

      for r_size in range(m + 1):
          for I_sub in itertools.combinations(range(m), r_size):
              # Calculate sum (using a tuple so it's hashable)
              current_sum = [0] * self.n
              for idx in I_sub:
                  vec = R_vecs[R_edges[idx]]
                  for k in range(self.n):
                      current_sum[k] += vec[k]

              # Store the first one we find for this sum
              sum_tuple = tuple(current_sum)
              if sum_tuple not in red_sum_map:
                  red_sum_map[sum_tuple] = I_sub

      # 3. Search for a matching Blue subset sum
      for b_size in range(n + 1):
          for J_sub in itertools.combinations(range(n), b_size):
              # Optimization: Skip the (0,0) and (m,n) case as per your logic
              # but we need to check the map for the blue sum

              blue_sum = [0] * self.n
              for idx in J_sub:
                  vec = B_vecs[B_edges[idx]]
                  for k in range(self.n):
                      blue_sum[k] += vec[k]

              blue_sum_tuple = tuple(blue_sum)

              if blue_sum_tuple in red_sum_map:
                  I_sub = red_sum_map[blue_sum_tuple]

                  # Check that we haven't just found the empty set or the full set
                  if (len(I_sub) == 0 and len(J_sub) == 0) or \
                     (len(I_sub) == m and len(J_sub) == n):
                      continue

                  E1 = set(R_edges[i] for i in I_sub) | set(B_edges[j] for j in J_sub)
                  E2 = set(all_incoming) - E1
                  return E1, E2

      return {}, {}



  def check_blue_range(self, args):
    """Worker function to check a specific range of blue sizes."""
    blue_indices, B_edges, B_vecs, red_sum_map, vector_len, m, n = args

    # blue_indices is a list of b_sizes to check, e.g., [1, 2]
    for b_size in blue_indices:
      for J_sub in itertools.combinations(range(n), b_size):
        blue_sum = [0] * vector_len
        for idx in J_sub:
          vec = B_vecs[B_edges[idx]]
          for k in range(vector_len):
            blue_sum[k] += vec[k]

        blue_sum_tuple = tuple(blue_sum)
        if blue_sum_tuple in red_sum_map:
          I_sub = red_sum_map[blue_sum_tuple]
          # Filter out the (0,0) and (m,n) forbidden cases
          if not ((len(I_sub) == 0 and len(J_sub) == 0) or
              (len(I_sub) == m and len(J_sub) == n)):
            return (I_sub, J_sub)
    return None

  def matching_partition_parallel(self, v):
    R_edges, B_edges, R_vecs, B_vecs = self.get_source_profile_vectors(v)
    m, n = len(R_edges), len(B_edges)

    # 1. Build the Red Sum Map (The "Meet-in-the-Middle" baseline)
    red_sum_map = {}
    for r_size in range(m + 1):
      for I_sub in itertools.combinations(range(m), r_size):
        current_sum = tuple(sum(x) for x in zip(*(R_vecs[R_edges[i]] for i in I_sub))) if I_sub else tuple([0]*self.n)
        if current_sum not in red_sum_map:
          red_sum_map[current_sum] = I_sub

    # 2. Prepare work for the Pool
    num_cores = cpu_count()
    # Split blue_sizes (0 to n) into chunks for each core
    blue_sizes = list(range(n + 1))
    chunk_size = max(1, len(blue_sizes) // num_cores)
    chunks = [blue_sizes[i:i + chunk_size] for i in range(0, len(blue_sizes), chunk_size)]

    worker_args = [
      (chunk, B_edges, B_vecs, red_sum_map, self.n, m, n)
      for chunk in chunks
    ]

    # 3. Execute in Parallel
    with Pool(processes=num_cores) as pool:
      # imap_unordered is slightly faster as it returns results as soon as they are ready
      for result in pool.imap_unordered(self.check_blue_range, worker_args):
        if result:
          I_sub, J_sub = result
          pool.terminate() # Stop other workers once we find a match
          E1 = set(R_edges[i] for i in I_sub) | set(B_edges[j] for j in J_sub)
          return E1, (set(R_edges) | set(B_edges)) - E1

    return {}, {}

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


  def matching_partition_claude(self, v):
    '''Get source profile vectors then find a proper partition of
    incoming edges into E1, E2 such that the sum of red vectors in E1
    equals the sum of blue vectors in E1.

    Optimizations vs original:
      - numpy for fast vector arithmetic
      - Precompute ALL blue subset sums into a hash dict (O(1) lookup)
      - Iterate only red subsets; each red_sum is a single dict lookup
      - Skip trivial (empty / full-set) partitions via set arithmetic
    '''
    R_edges, B_edges, R_vecs, B_vecs = self.get_source_profile_vectors(v)
    all_incoming = list(set(R_edges) | set(B_edges))
    if len(all_incoming) < 2:
      return None

    # --- Shortcut: isolate any edge whose vector is all zeros ---
    for edge in all_incoming:
      vec = R_vecs.get(edge) or B_vecs.get(edge)
      if not any(vec):            # faster than all(x==0 for x in vec)
        E1 = {edge}
        return E1, set(all_incoming) - E1

    m, n = len(R_edges), len(B_edges)

    # --- Convert to numpy arrays once ---
    R_arrs = [np.array(R_vecs[e], dtype=np.int64) for e in R_edges]
    B_arrs = [np.array(B_vecs[e], dtype=np.int64) for e in B_edges]
    zero   = np.zeros(self.n, dtype=np.int64)

    # --- Precompute every blue subset sum (2ⁿ entries) ---
    # Maps tuple(blue_sum) -> (r_indices, b_indices) for the first hit found.
    # We record b_indices here; r_indices gets filled in the red loop.
    blue_sums: dict[tuple, tuple] = {}
    for b_size in range(n + 1):
      for J_sub in itertools.combinations(range(n), b_size):
        b_sum = sum((B_arrs[j] for j in J_sub), zero.copy())
        key   = tuple(b_sum)
        if key not in blue_sums:      # keep first (smallest) match
          blue_sums[key] = J_sub

    # --- Iterate red subsets; O(1) lookup into blue_sums ---
    for r_size in range(m + 1):
      for I_sub in itertools.combinations(range(m), r_size):
        r_sum = sum((R_arrs[i] for i in I_sub), zero.copy())
        J_sub = blue_sums.get(tuple(r_sum))
        if J_sub is None:
          continue

        E1 = (
          {R_edges[i] for i in I_sub} |
          {B_edges[j] for j in J_sub}
        )
        # Must be a *proper* subset (not empty, not everything)
        if E1 and E1 != set(all_incoming):
          return E1, set(all_incoming) - E1

    return {}, {}
