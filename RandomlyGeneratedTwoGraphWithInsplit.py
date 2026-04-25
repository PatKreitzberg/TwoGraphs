from RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph

class RandomlyGeneratedTwoGraphWithInsplit(RandomlyGeneratedTwoGraph):
  def __init__(self, n, z):
    super().__init__(n,v)

  def force_insplit(self):
    pass


  def force_insplit_partition(self):
    '''
    Want a vertex v such that there is a vertex u with R*B[u][v] = R*B[u][v] = 0
    and
    there are two vertices a,b such that (R+B)[a][v] = (R+B)[b][v] = 0

    Then we can add edges:

    u->a->v
    u->b->v
    to get commuting square
    '''

    for row in range(self.n): # source vertex
      for col in range(self.n): # range vertex
        assert len(self.RB_paths_matrix[row][col]) == len(self.BR_paths_matrix[row][col])
    print("PASSED ASSERT")

    u,a,b = None,None,None
    for v in range(self.n):
      u = self.v_has_u(v)
      if u is None:
        continue

      a,b = self.v_has_ab(v)
      if a is None and b is None:
        continue

      # v is insplittable!
      self.insplit_v = v
      self.can_insplit = True
      break

    if (u is None) or (a is None) or (b is None):
      return None,None

    #add edges
    # 1. u -> a degree 1
    # 2. u -> b degree 2
    # 3. b -> v degree 1
    # 4. a -> v degree 1
    # E1 = {b->v, a->v}

    # Red edges
    ua_edge = PathMatrix.create_edge(u, a, self.R_degree, 0)
    bv_edge = PathMatrix.create_edge(b, self.insplit_v, self.R_degree, 0)
    self.R_path_matrix.add_edge(ua_edge)
    self.R_path_matrix.add_edge(bv_edge)

    # Blue edges
    ub_edge = PathMatrix.create_edge(u, b, self.B_degree, 0)
    av_edge = PathMatrix.create_edge(a, self.insplit_v, self.B_degree, 0)
    self.B_path_matrix.add_edge(ub_edge)
    self.B_path_matrix.add_edge(av_edge)

    self.RB_paths_matrix = self.R_path_matrix*self.B_path_matrix
    self.BR_paths_matrix = self.B_path_matrix*self.R_path_matrix

    for row in range(self.n): # source vertex
      for col in range(self.n): # range vertex
        assert len(self.RB_paths_matrix[row][col]) == len(self.BR_paths_matrix[row][col])
    print("PASSED ASSERT")



  def v_has_u(self, v):
    for u in range(self.n):
      if len(self.RB_paths_matrix[u][v]) == 0:
        return u
    return None

  def v_has_ab(self, v):
    ab = set()
    for i in range(self.n):
      if len(self.R_path_matrix[i][v]) + len(self.B_path_matrix[i][v]) == 0:
        ab.add(i)
        if len(ab) == 2:
          return ab
    return None,None


  def get_commuting_squares_for_partition(self, E_red_edges, E_blue_edges, R, B):
    red_blue_paths = set()
    for red_edge in E_red_edges: # Red edges
      # it is the range edge so we need to look at the other matrix to get source edge
      for source_vertex in range(B.n):
        for blue_edge in B[source_vertex][red_edge.s]:
          red_blue_paths.add((red_edge, blue_edge))

    blue_red_paths = set()
    for blue_edge in E_blue_edges:
      for source_vertex in range(R.n):
        for red_edge in R[source_vertex][blue_edge.s]:
          blue_red_paths.add((blue_edge, red_edge))
    commuting_squares = set()

    assert len(red_blue_paths) == len(blue_red_paths)

    for (r1,s1), (r2,s2) in zip(red_blue_paths, blue_red_paths):
      print("sources:", s1.s, s2.s)
      print("ranges:", r1.r, r2.r)
      commuting_squares.add(CommutingSquare(r1,s1,r2,s2))

    return commuting_squares

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
