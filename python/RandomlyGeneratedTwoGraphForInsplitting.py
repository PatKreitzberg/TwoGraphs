import numpy as np
import itertools
import collections
from collections import defaultdict as dd
import multiprocessing
from functools import partial
import time,datetime

from python.CommutingSquare import CommutingSquare
from python.Edge import Edge
from python.RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph

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


class RandomlyGeneratedTwoGraphForInsplitting(RandomlyGeneratedTwoGraph):
  def __init__(self, n, z, R=None, B=None, symmetric=False):
    super().__init__(n,z, R=R, B=B, symmetric=symmetric)

    self.time_limit_seconds = datetime.timedelta(seconds=10)
    # Find vertex which can be partitioned at
    v,E1,E2 = self.find_insplit_vertex_opt()

    self.is_legit = True
    if v is None:
      self.is_legit = False
      return

    assert v is not None
    assert len(E1) > 0
    assert len(E2) > 0

    if type(self) is RandomlyGeneratedTwoGraphForInsplitting:
      # This means we are not inheriting
      commuting_squares = self.get_commuting_squares(v, E1,E2)
      assert commuting_squares is not None

      n_commuting_squares = len(commuting_squares)

      n_rb_paths = sum([sum([len(col) for col in row]) for row in self.BR_paths_matrix])
      assert n_commuting_squares == n_rb_paths

      self.v = v
      self.E1 = E1
      self.E2 = E2

      self.vertices = [i for i in range(n)]
      self.edges = {edge for edge in self.R_path_matrix.edges + self.B_path_matrix.edges}
      self.commuting_squares = commuting_squares
      if len(commuting_squares) == 0:
        self.is_legit = False
        return

      self.calculate_boundary_matrices()


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

        commuting_squares = self.commuting_squares_for_vertices(row, col, commuting_squares)

    # Commuting squares for the insplit vertex
    E1_RB_paths,E1_BR_paths = set(),set()
    E2_RB_paths,E2_BR_paths = set(),set()

    for row in range(self.n):
      E1_RB_paths, E2_RB_paths = self.partition_paths(E1, E2, self.RB_paths_matrix[row][v])
      E1_BR_paths, E2_BR_paths = self.partition_paths(E1, E2, self.BR_paths_matrix[row][v])

      for E_RB_paths,E_BR_paths in [(E1_RB_paths, E1_BR_paths), (E2_RB_paths, E2_BR_paths)]:
        for RB_path, BR_path in zip(E_RB_paths,E_BR_paths):
          commuting_squares = self.commuting_square_for_vertices(RB_path, BR_path, commuting_squares)

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



  def matching_partition_opt(self, v):
    time_start = datetime.datetime.now()

    R_edges, B_edges, R_vecs, B_vecs = self.get_source_profile_vectors(v)
    all_incoming = list(set(R_edges) | set(B_edges))

    if len(all_incoming) < 2:
      return set(), set()

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
        time_elapsed = datetime.datetime.now() - time_start
        if time_elapsed > self.time_limit_seconds:
          print('taking too long')
          return set(), set()

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
        time_elapsed = datetime.datetime.now() - time_start
        if time_elapsed > self.time_limit_seconds:
          print('taking too long')
          return set(), set()

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


  def find_insplit_vertex_opt(self):
    v_to_graph_degree = {v:0 for v in range(self.n)}
    for se in range(self.n):
      for re in range(self.n):
        v_to_graph_degree[re] += self.R_path_matrix.adj_matrix[se][re]
        v_to_graph_degree[re] += self.B_path_matrix.adj_matrix[se][re]

    #v_sorted_by_degree = sorted( [(int(d),v) for v,d in v_to_graph_degree.items()] )
    v_sorted_by_degree =  [(int(d),v) for v,d in v_to_graph_degree.items()]

    for _,v in v_sorted_by_degree:
      E1, E2 =  self.matching_partition_opt(v)
      if len(E1)>0 and len(E2)>0:
        return v, E1, E2
    return None,{},{}

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
