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
          if (len(self.path_matrix[Ar][ell]) > 0) and (len(other.path_matrix[ell][Bc]) > 0):
            for A_edge in self.path_matrix[Ar][ell]:
              for B_edge in other.path_matrix[ell][Bc]:
                # Edges are of the form
                # (color, source vertex, range vertex, edge key)
                # Edge key is just if there are m-many 'r' edges between 0,1 we get
                # ('r',0,1,0), ('r',0,1,1), ..., ('r',0,1, m)
                print("A_edge", A_edge)
                print("B_edge", B_edge)


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
        [
          (color, (r,c), edge_key) for edge_key in range(A[r][c])
        ]
        for c in range(self.n)
    ]
      for r in range(self.n)
    ]
    return path_matrix


   # path_matrix =  [
   #   [
   #     ((color, r,c), A[r][c]) for c in range(self.n)
   # ]
   #   for r in range(self.n)
   # ]
   # return path_matrix
