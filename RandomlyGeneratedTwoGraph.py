from TwoGraph import TwoGraph
from PathMatrix import PathMatrix
from CommutingSquare import CommutingSquare

class RandomlyGeneratedTwoGraph(TwoGraph):
  def __init__(self, n, z):
    # n is number of vertices
    # z is an upper limit on number of edges in graph between any two (maybe non-distinct) vertices

    R_degree = 1
    B_degree = 2
    R,B = self.generate_adjacency_matrices(n, z)

    R_path_matrix = PathMatrix(R, R_degree)
    B_path_matrix = PathMatrix(B, B_degree)

    edges = R_path_matrix.edges + B_path_matrix.edges
    commuting_squares = self.get_commuting_squares(R_path_matrix, B_path_matrix)
    print("Okay creating TwoGraph now")
    load_from = {
      'vertices':[i for i in range(n)],   # Vertices are just basic integers 0,1,...,n-1
      'edge_label_to_edge':{edge.label:edge for edge in edges},
      'commuting_squares':commuting_squares
    }
    super().__init__(load_from)

  def generate_adjacency_matrices(self, n, z):
    ''''
    Generates adjacency matrices R and B
    R is just for red edges, B for blue edges

    CONSTRAINTS:
    RB = BR so that for every length two red-blue path there is an equivalent blue-red path
    Elements of R, B are nonnegative integers between 0 and z, roughly. For the first generated matrix this is
    always true; for the scond matrix it is harder to guarantee because it must be tailored to commute with
    the first matrix
    '''

    R = [[1,0],
           [0,1]]
    B = [[0,2],  # two eges from 0 to 1
           [0,0]]  # A[row][column]

    return R,B

  def get_commuting_squares(self, R, B):
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
    RB_paths_matrix = R*B
    BR_paths_matrix = B*R
    commuting_squares = []

    n = R.n
    for i in range(n):
      for j in range(n):
       assert len(RB_paths_matrix[i][j]) == len(BR_paths_matrix[i][j])

       for RB_path, BR_path  in zip(RB_paths_matrix[i][j], BR_paths_matrix[i][j]):
         s1,r1 = RB_path
         s2,r2 = BR_path
         commuting_squares.append(CommutingSquare(r1,s1, r2, s2))
    return commuting_squares
