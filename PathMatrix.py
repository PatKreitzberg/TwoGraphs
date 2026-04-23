class PathMatrix:
  def __init__(self, A, color):
    '''
    A: Adjacency matrix for DIRECTED GRAPH
    Returns a matrix from A where each entry is (color, (i,j), A_ij)
    This allows us to get the path
    '''
    self.n = len(A)
    self.path_matrix =  self.get_path_matrix(A, color)


  def get_path_matrix(self, A, color):
    '''
    A: Adjacency matrix for DIRECTED GRAPH
    Returns a matrix from A where each entry is (color, (i,j), A_ij)
    This allows us to get the path


    (color, (i,j), number of edges from i -> j)
    '''
    return [
      [
        (color, (r,c), A[r][c]) for r in range(n)
      ]
      for c in range(n)
    ]

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
          (Acolor, (Ai,Aj), Ak) = self.path_matrix[Ar][ell]  # Move across row in A
          (Bcolor, (Bi,Bj), Bk) = other.path_matrix[ell][Bc] # Move down column in B

          # (Ai, Aj) is an edge from v_i to v_j in the A matrix
          # Ak is number of such edges (Ak = 0 if no edges exist)

          # Need to add new path for Ak*Bk because there are Ak many
          # edges from Ai to Aj and Bk many edges from Bi to Bj. So
          # there are Ak*Bk many paths from Ai to Bj.

          # if there are no paths from vi to vj in A then Ak = 0 so there are none appended to res[Ar][Bc]

          color_of_paths = (Acolor, Bcolor)
          for pi in range(Ak*Bk):
            path = (Ai, Aj, Bi, Bj)
            index = pi
            res[Ar][Bc].append( (color_of_paths, path, index) )
    return res


R = [[1,0],
       [0,1]]
B = [[0,2],  # two eges from 0 to 1
       [0,0]]  # A[row][column]

n = len(R)


#print(R)
print()
pR = PathMatrix(R, 'r')
pB = PathMatrix(B, 'b')
res = pR*pB

st = ''
for r in res:
  for i in range(len(r)):
    for p in r[i]:
      if p[-1] > 0:
        print(p)


print('pR')
print(pR)
print('pB')
print(pB)
