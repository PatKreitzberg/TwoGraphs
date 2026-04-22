from sympy import Matrix,pprint
from sympy.matrices.normalforms import smith_normal_form, hermite_normal_form
from fractions import Fraction
import math

class BoundaryFunctionMatrix:
  def __init__(self, graph, r, domain_items, range_items, domain_item_to_index, range_item_to_index, calc_ker=False, calc_img=False):
    self.r = r
    self.matrix = self.build_matrix(graph, domain_items, range_items, domain_item_to_index, range_item_to_index)

    self.img = None
    self.img_str = None

    self.ker = None
    self.ker_str = None

    if calc_img:
      self.calc_img()
      self.img_str = self.partial_span_as_str(self.img, range_items)

    if calc_ker:
      self.calc_ker()
      self.ker_str_items = self.partial_span_as_str(self.ker, domain_items)
      self.ker_str = f"Z^{len(self.ker)}"

  def calc_ker(self):
    if self.ker is None:
      self.ker = self.kernel_basis(self.matrix)
    return self.ker

  def calc_img(self):
    if self.img is None:
      self.img = self.image_basis(self.matrix)
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

  def build_matrix(self, graph, domain_items, range_items, domain_item_to_index, range_item_to_index):
    # if r = 1
    # domain items are edges
    # range items are vertices

    # if r = 2
    # domain items are commuting squares
    # range items are edges

    n_range_items = len(range_items)
    n_domain_items = len(domain_items)
    matrix = [[0]*n_domain_items for _ in range(n_range_items)] # matrix of all zeros

    for i in range(1, self.r+1):
      for domain_item in domain_items:
        domain_item_index = domain_item_to_index[domain_item]
        for ell in [0,1]:
          range_item_index = range_item_to_index[domain_item.F(i,ell)]
          #print(domain_items[domain_item_index], range_items[range_item_index], f"F_{i}^{ell} = {range_items[range_item_to_index[domain_item.F(i,ell)]]}")
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
            span_vec_str += str(items[item_index])
          elif span_vec[item_index] == -1:
            span_vec_str += '-' + str(items[item_index])
          else:
            span_vec_str += str(span_vec[item_index]) + str(items[item_index])
          span_vec_str += ' + '

      if span_vec_str.endswith(' + '):
        span_vec_str = span_vec_str[:-2]
      #print(span_vec_str)
      out += span_vec_str + ', '
    if out == "":
      return "{0}"
    if out.endswith(', '):
      out = out[:-2].strip()
    return out
