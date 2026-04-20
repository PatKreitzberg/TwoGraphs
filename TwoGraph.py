import re
from collections import defaultdict as dd
from fractions import Fraction
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

from EdgeAndCommutingSquare import Edge, CommutingSquare
from MatrixMath import PartialOneMatrix

class TwoGraph:
  def __init__(self, path):
    self.edge_label_to_edge = {}
    self.edge_to_degree = {}
    self.degree_to_edges = dd(list)
    self.commuting_squares = []
    self.edge_to_commuting_squares = dd(set)
    self.vertices = set()
    self.parse(path)
    self.edges = self.edge_label_to_edge.keys()
    self.vertex_to_index = {v:i for i,v in enumerate(self.vertices)}
    self.edge_to_index = {e:i for i,e in enumerate(self.edges)}

    # boundary function matrices
    self.partial_one = PartialOneMatrix(self)
    print("Partial one matrix:", self.partial_one.matrix)

  # def image_basis(self, A: list[list[int]]) -> list[list[int]]:
  #     """
  #     Return an integer basis for the image (column space) of matrix A.

  #     Parameters
  #     ----------
  #     A : 2D list of ints, shape (m, n).

  #     Returns
  #     -------
  #     List of column vectors (each a list of ints) forming a ℤ-basis for
  #     im(A) ⊆ ℤ^m.  Returns an empty list when A is the zero matrix.
  #     """
  #     M = Matrix(A)
  #     # HNF(M) is upper-triangular with zero columns dropped,
  #     # so its columns are exactly a ℤ-basis for im(A).
  #     H = hermite_normal_form(M)
  #     basis = []
  #     for j in range(H.cols):
  #         col = [int(H[i, j]) for i in range(H.rows)]
  #         if any(x != 0 for x in col):
  #             basis.append(col)
  #     return basis


  # def kernel_basis(self, A: list[list[int]]) -> list[list[int]]:
  #     """
  #     Return an integer basis for the kernel (null space) of matrix A.

  #     Parameters
  #     ----------
  #     A : 2D list of ints, shape (m, n).

  #     Returns
  #     -------
  #     List of column vectors (each a list of ints) forming a ℤ-basis for
  #     ker(A) ⊆ ℤ^n.  Returns an empty list when the kernel is trivial.
  #     """
  #     M = Matrix(A)
  #     rational_ns = M.nullspace()   # exact rational arithmetic

  #     basis = []
  #     for v in rational_ns:
  #         # Scale by LCM of denominators to clear fractions.
  #         fracs = [Fraction(str(entry)) for entry in v]
  #         denom_lcm = 1
  #         for f in fracs:
  #             denom_lcm = denom_lcm * f.denominator // math.gcd(denom_lcm, f.denominator)
  #         basis.append([int(f * denom_lcm) for f in fracs])
  #     return basis



  def print_commuting_squares(self):
    # Print commuting squares
      for cs in self.commuting_squares:
        print("Commuting square:", cs)
        for i in [1,2]:
          for ell in [0,1]:
            print(f"F_{i}^{ell} = {cs.F(i,ell).edge_label}" )

  def parse(self, file_path):
    current_section = None
    with open(file_path, 'r') as f:
      for line in f:
        line = line.strip()

        # Skip empty lines
        if not line:
          continue

        # Detect section headers
        if line.startswith('#'):
          header = line.lower()
          if 'edges' in header:
            current_section = 'edges'
          elif 'degrees' in header:
            current_section = 'degrees'
          elif 'commuting squares' in header:
            current_section = 'commuting_squares'
          continue

        # Parse based on the current section
        if current_section == 'edges':
          # Format: <label> <v1> <v2>
          parts = line.split()
          if len(parts) == 3:
            e = parts[0]
            s = parts[1]
            r = parts[2]
            edge = Edge(e,s,r)
            self.vertices.add(s)
            self.vertices.add(r)
            assert not (e in self.edge_label_to_edge.keys())
            self.edge_label_to_edge[e] = edge

        elif current_section == 'degrees':
          parts = line.split()
          if len(parts) > 1:
            degree_index = int(parts[0])
            for e in parts[1:]:
              self.degree_to_edges[degree_index].append(self.edge_label_to_edge[e])
              self.edge_label_to_edge[e].degree_index = degree_index

        elif current_section == 'commuting_squares':
         # Format: label_a label_b = label_c label_d
         if '=' in line:
           left_side, right_side = line.split('=')
           lhs = tuple(left_side.strip().split())
           rhs = tuple(right_side.strip().split())

           left_edge_1 = self.edge_label_to_edge[lhs[0]]
           left_edge_2 = self.edge_label_to_edge[lhs[1]]
           right_edge_1 = self.edge_label_to_edge[rhs[0]]
           right_edge_2 = self.edge_label_to_edge[rhs[1]]

           # ab ~ cd
           commuting_square = CommutingSquare(left_edge_1, left_edge_2, right_edge_1, right_edge_2)
           self.commuting_squares.append(commuting_square)
           for e in [left_edge_1, left_edge_2, right_edge_1, right_edge_2]:
             self.edge_to_commuting_squares[e].add(commuting_square)
