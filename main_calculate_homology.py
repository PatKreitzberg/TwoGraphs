from sage.all import vector, matrix, ZZ, ChainComplex
from python.TwoGraph import TwoGraph

def cohomology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  cochain = ChainComplex({1: d1, 2: d2}, degree=-1,base_ring=ZZ).dual()
  C0 = cochain.homology(0)
  C1 = cochain.homology(1)
  C2 = cochain.homology(2)
  return C0,C1,C2

def homology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1, base_ring=ZZ)
  H0 = C.homology(0)
  H1 = C.homology(1)
  H2 = C.homology(2)
  return H0,H1,H2


def check_cohomology(gfile):
  g = TwoGraph(gfile)
  C0,C1,C2 = cohomology(g)
  print("Standard cohomology")
  print('C0', C0)
  print('C1', C1)
  print('C2', C2)
  delta_0 = matrix(ZZ, g.d_1.matrix).transpose()
  delta_1 = matrix(ZZ, g.d_2.matrix).transpose()
  cochain = ChainComplex({0: delta_0, 1: delta_1}, degree=1,base_ring=ZZ)
  C0 = cochain.homology(0)
  C1 = cochain.homology(1)
  C2 = cochain.homology(2)
  print("New cohomology")
  print('C0', C0)
  print('C1', C1)
  print('C2', C2)
