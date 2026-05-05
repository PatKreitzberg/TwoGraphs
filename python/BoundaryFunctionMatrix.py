from sympy import Matrix,pprint
from sympy.matrices.normalforms import smith_normal_form, hermite_normal_form
from fractions import Fraction
import math
from sage.matrix.constructor import matrix

from python.CommutingSquare import CommutingSquare

class BoundaryFunctionMatrix:
  def __init__(self,
               graph,
               r,
               domain_items,
               range_items,
               domain_item_to_index,
               range_item_to_index,
               calc_ker=False,
               calc_img=False):

    self.r = r
    self.domain_items = domain_items
    self.range_items = range_items
    self.domain_item_to_index = domain_item_to_index
    self.range_item_to_index = range_item_to_index

    self.matrix = self.build_matrix(graph)

    self.img = None
    self.img_str = None

    self.ker = None
    self.ker_str = None

    if calc_img:
      self.calc_img()

    if calc_ker:
      self.calc_ker()

  def calc_ker(self):
    if self.ker is None:
      A = matrix(self.matrix)
      self.ker =  [list(v) for v in A.right_kernel().basis()]
    self.ker_str_items = self.partial_span_as_str(self.ker, self.domain_items)
    self.ker_str = f"Z^{len(self.ker)}"
    return self.ker

  def calc_img(self):
    if self.img is None:
      A = matrix(self.matrix)
      self.img = [list(v) for v in A.column_space().basis()]
      self.img_str_items = self.partial_span_as_str(self.img, self.range_items)
    return self.img

  def __str__(self):
    out = ''
    for row in self.matrix:
      out += str(row) + '\n'
    return out

  def increment(self, i, ell):
    if (i + ell)%2 == 0:
      return 1
    return -1

  def build_matrix(self, graph):
    # if r = 1
    # domain items are edges
    # range items are vertices

    # if r = 2
    # domain items are commuting squares
    # range items are edges

    n_range_items = len(self.range_items)
    n_domain_items = len(self.domain_items)
    self.matrix_dim = (n_range_items, n_domain_items)

    # one column for each domain item
    # one row for each range item
    matrix = [[0]*n_domain_items for _ in range(n_range_items)] # matrix of all zeros

    for i in range(1, self.r+1):
      for domain_item in self.domain_items:
        domain_item_index = self.domain_item_to_index[domain_item]

        for ell in [0,1]:
          range_item_index = self.range_item_to_index[domain_item.F(i,ell)]
          matrix[range_item_index][domain_item_index] += self.increment(i,ell)
    return matrix

  def image_basis(self, A: list[list[int]]) -> list[list[int]]:
      M = Matrix(A)
      # HNF(M) is upper-triangular with zero columns dropped,
      # so its columns are exactly a ℤ-basis for im(A).
      H = hermite_normal_form(M)
      basis = []
      for j in range(H.cols):
          col = [int(H[i, j]) for i in range(H.rows)]
          if any(x != 0 for x in col):
              basis.append(col)
      return basis


  def kernel_basis(self, A: list[list[int]]) -> list[list[int]]:
      M = Matrix(A)
      rational_ns = M.nullspace()   # exact rational arithmetic

      basis = []
      for v in rational_ns:
          # Scale by LCM of denominators to clear fractions.
          fracs = [Fraction(str(entry)) for entry in v]
          denom_lcm = 1
          for f in fracs:
              denom_lcm = denom_lcm * f.denominator // math.gcd(denom_lcm, f.denominator)
          basis.append([int(f * denom_lcm) for f in fracs])
      return basis

  def partial_span_as_str(self, part, items):
    out = ''
    for span_vec in part:
      span_vec_str = ''

      for item_index in range(len(span_vec)):
        if span_vec[item_index] != 0:
          if span_vec[item_index] == 1:
            span_vec_str +=  '$' + str(items[item_index]) + '$'

          elif span_vec[item_index] == -1:
            span_vec_str += '$-' + str(items[item_index]) + '$'

          else:
            span_vec_str +=   '$'+str(span_vec[item_index]) + str(items[item_index]) + '$'

          span_vec_str += ' + '
      if span_vec_str.endswith(' + '):
        span_vec_str = span_vec_str[:-2]
      out += span_vec_str + ', '
    if out == "":
      return "{0}"
    if out.endswith(', '):
      out = out[:-2].strip()
    return out


  def latex(self):
    n_rows = len(self.matrix) # range items
    n_cols  = len(self.matrix[0])     # domain items
    nl = '\n'
    start_of_matrix = '\\[' + nl
    start_of_matrix += '\\begin{blockarray}{r' + 'r'*(n_cols) + '}' + nl

    end_of_matrix = '\\end{block}\n\\end{blockarray}\n\\]\n'
    # actual matrix
    index_to_domain_item = {i:d for d,i in self.domain_item_to_index.items()}
    index_to_range_item = {i:d for d,i in self.range_item_to_index.items()}

    # Label the columns, either as edges for r=1 or commuting squares for r=2
    matrix_col_labels = '   & '
    for ci in range(n_cols):
      if type(self.domain_items[0]) is CommutingSquare:
        matrix_col_labels +=self.domain_items[self.domain_item_to_index[index_to_domain_item[ci]]].latex_label + ' & '
      else:
        matrix_col_labels +=str(self.domain_items[self.domain_item_to_index[index_to_domain_item[ci]]]) + ' & '

    matrix_col_labels = matrix_col_labels[:-2] + '\\\\' + nl
    start_of_matrix += matrix_col_labels
    start_of_matrix += '\\begin{block}{r(' + 'r'*(n_cols) + ')}' + nl
    matrix = ''
    for ri in range(n_rows):
      matrix += str(self.range_items[self.range_item_to_index[index_to_range_item[ri]]])
      matrix += ' & '
      for ci in range(n_cols):
        matrix += str(self.matrix[ri][ci]) + ' & '
      matrix = matrix[:-2] + '\\\\' + nl
    return start_of_matrix + matrix + end_of_matrix
