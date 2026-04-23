class PathMatrix:
  def __init__(self, A, color):
    '''
    A: Adjacency matrix for DIRECTED GRAPH
    Returns a matrix from A where each entry is (color, (i,j), A_ij)
    This allows us to get the path
    '''
    self.n = len(A)
    print("Got matrix", A)
    self.path_matrix =  self.get_path_matrix(A, color)

  def __getitem__(self, i):
    return self.path_matrix[i]

  def __iter__(self):
    return iter(self.path_matrix)

  def __str__(self):
    out = 'Path matrix:\n'
    for row in self.path_matrix:
      for col in row:
        out += str(col) + ' '
      out += '\n'
    return out

  def __mul__(self, other):
    '''
    If we have path matrices A, B and this one is A then __mul__ is
    called when user has A*B. If want to allow for B*A then do
    __rmul__ but if B has a __mul__ then that is called. this is only
    the backup

     Calculates  A * B
    Inputs:
     A is self.path_matrix
     B is other.path_matrix

    Output:
     res is n x n matrix where each entry (i,j) is a list of color A - color B paths from vi to vj
       if there are no paths then the entry (i,j) is just an empty list
    '''
    assert len(other.path_matrix) == self.n

    res = [[[] for i in range(self.n)] for j in range(self.n)]

    for Ar in range(self.n): #for row in A
      for Bc in range(self.n):
        paths = []
        for ell in range(self.n):
          A_edge, Ak = self.path_matrix[Ar][ell]  # Move across row in A
          B_edge, Bk = other.path_matrix[ell][Bc] # Move down column in B

          if (Ak*Bk) > 0:
            print("A edge", A_edge, Ar, ell)
            print("B edge", B_edge, ell, Bc)

          # edges are of the form (Color, i, j) which is a colored edge from vertex i to vertex j
          # Ak is number of such edges (Ak = 0 if no edges exist)

          # Need to add new path for Ak*Bk because there are Ak many
          # edges from Ai to Aj and Bk many edges from Bi to Bj. So
          # there are Ak*Bk many paths from Ai to Bj.

          # if there are no paths from vi to vj in A then Ak = 0 so there are none appended to res[Ar][Bc]

          path = (A_edge, B_edge)
          res[Ar][Bc] += [
            (path, path_index) for path_index in range(Ak*Bk)
          ]

    return res


  def get_path_matrix(self, A, color):
    '''
    A: Adjacency matrix for DIRECTED GRAPH
    Returns a matrix from A where each entry is (color, (i,j), A_ij)
    This allows us to get the path


    (color, (i,j), number of edges from i -> j)
    '''
    path_matrix =  [
      [
        ((color, r,c), A[r][c]) for c in range(self.n)
    ]
      for r in range(self.n)
    ]
    return path_matrix
