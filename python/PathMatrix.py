from python.Edge import Edge

class PathMatrix:
  def __init__(self, A, degree):
    '''
    A: Adjacency matrix for DIRECTED GRAPH
    Returns a matrix from A with A(i,j) is a list of edges from vertex i to vertex j
    Vertices are integers {0,1,...,n-1}
    Edge labels are:
      E(deg=<degree>, (source vertex, range vertex), #=<edge key>)
      edge key: there may be multiple edges between source and range of this degree so we separate using this number
    '''
    self.n = len(A)
    self.degree = degree
    self.adj_matrix = A
    self.path_matrix =  self.get_path_matrix(A, degree)
    self.edges = self.get_edges()
    self.size = sum([sum(r) for r in self.adj_matrix])

  def add_edge(self, edge):
    self.path_matrix[edge.s][edge.r].append(edge)
    self.adj_matrix[edge.s][edge.r] += 1
    self.edges.append(edge)

  def __getitem__(self, i):
    return self.path_matrix[i]

  def __iter__(self):
    return iter(self.path_matrix)

  def __str__(self):
    out = 'Path matrix:\n'
    for row in self.path_matrix:
      for col in row:
        out += str([e.label for e in col]) + ' '
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
     res is n x n matrix where each entry (i,j) is a list of degree A - degree B paths from vi to vj
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
                # (degree, source vertex, range vertex, edge key)
                # Edge key is just if there are m-many 'r' edges between 0,1 we get
                # ('r',0,1,0), ('r',0,1,1), ..., ('r',0,1, m)
                res[Ar][Bc].append((A_edge, B_edge))

    return res


  def get_path_matrix_small_n(self, A, degree):
    self.edge_labels_superset = set()
    if degree == 1: # red edges
      self.edge_labels_superset = [s for s in 'abcdefghij']
    elif degree == 2:
      self.edge_labels_superset = [s for s in 'mnpqrstwxyz']
    else:
      assert False

    # Reverse so a is poped first then b...
    self.edge_labels_superset = self.edge_labels_superset[::-1]
    path_matrix =  [
      [
        [
          self.create_edge_small_n(row, col, self.degree, edge_key) for edge_key in range(A[row][col])
        ]
        for col in range(self.n)
      ]
      for row in range(self.n)
    ]
    return path_matrix

  def create_edge_small_n(self, source_vertex, range_vertex, degree, edge_key):
    edge_label = self.edge_labels_superset.pop()
    return Edge(edge_label, source_vertex, range_vertex, degree=degree)


  def get_path_matrix(self, A, degree):
    '''
    A: Adjacency matrix for DIRECTED GRAPH
    Returns a matrix from A where each entry is (degree, (i,j), A_ij)
    This allows us to get the path
    Edges will be Edge items with
    label = "E(row,col)edge_key"
    s = row
    r = col
    degree = degree
    '''

    total_edges = 0
    for row in A:
      total_edges += sum(row)
    if total_edges <= 10:
      return self.get_path_matrix_small_n(A, degree)

    # Edge of the form
    path_matrix =  [
      [
        [
          PathMatrix.create_edge(row, col, self.degree, edge_key) for edge_key in range(A[row][col])
        ]
        for col in range(self.n)
    ]
      for row in range(self.n)
    ]
    return path_matrix

  def get_edges(self):
    edges = []
    for r in self.path_matrix:
      for eset in r:
        edges += eset
    return edges

  @staticmethod
  def create_edge(source_vertex, range_vertex, degree, edge_key):
    edge_label = 'E(deg=' + str(degree) + ',(' + str(source_vertex) + ',' + str(range_vertex) + ')#=' + str(edge_key)
    return Edge(edge_label, source_vertex, range_vertex, degree=degree)
