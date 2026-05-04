import sys, random
import numpy as np

from RandomlyGeneratedTwoGraphForInsplitting import RandomlyGeneratedTwoGraphForInsplitting
from InsplitTwoGraph import InsplitTwoGraph
from RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph
from TwoGraph import TwoGraph

from sage.all import matrix, ZZ, ChainComplex

from pprint import pprint


def sage_cohomology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  cochain = ChainComplex({1: d1, 2: d2}, degree=-1,base_ring=ZZ).dual()
  C0 = cochain.homology(0)
  C1 = cochain.homology(1)
  C2 = cochain.homology(2)
  return C0,C1,C2

def my_cohomology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  cochain = ChainComplex({1: d1, 2: d2}, degree=-1,base_ring=ZZ)


def inspect_cohomology(g):
  delta1 = matrix(ZZ, g.d_1.matrix).transpose()
  delta2 = matrix(ZZ, g.d_2.matrix).transpose()

  return


def sage_homology(g, generators=False):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1,base_ring=ZZ)
  #print('c.var')
  #print(vars(C))
  #print()
  #print('dir(C)')
  #for x in dir(C):
  #  print(x)
  #print()
  #print('differential', C.differential())

  H0 = C.homology(0, generators=generators)
  H1 = C.homology(1, generators=generators)
  H2 = C.homology(2, generators=generators)
  return H0,H1,H2

def print_generators(title, H0,H1,H2):
  print(title)
  print('H0:',H0)
  print()
  print('H1:',H1)
  print()
  print('H2:',H2)
  print()

def print_snf(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  D1,U1,V1 = d1.smith_form(integral=True)

  print("Smith normal forms:")
  print("Partial_1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print('d1 matrix')
  print(d1)
  print("D")
  print(D1)
  print("U1*d1*V1")
  print(U1*d1*V1)
  assert D1 == U1*d1*V1
  assert d1 == (U1.inverse())*D1*(V1.inverse())

  d2 = matrix(ZZ, g.d_2.matrix)
  D2,U2,V2 = d2.smith_form(integral=True)
  print("Partial_2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("d2")
  print(d2)
  print(D2)
  print("U2*d2*V2")
  print(U2*d2*V2)


  assert D2 == U2*d2*V2
  assert d2 == (U2.inverse())*D2*(V2.inverse())



def print_homology(g, gi, cohomology=False):
  print("Original graph~~~~~~~~~")
  if cohomology:
    print("\tCohomology")
    H0,H1,H2 = sage_cohomology(g)
    print()
    inH0,inH1,inH2 = sage_cohomology(gi)
  else:
    print("\tHomology")
    H0,H1,H2 = sage_homology(g)
    print()
    inH0,inH1,inH2 = sage_homology(gi)

  print("\tH0", H0)
  print("\tH1", H1)
  print("\tH2", H2)
  print("Insplit graph~~~~~~~~~")
  print("\tH0", inH0)
  print("\tH1", inH1)
  print("\tH2", inH2)


def snf(g, gi):
  g = TwoGraph(g)
  gi = TwoGraph(gi)
  print("Original graph~~~~~~~~~")
  print_snf(g)
  print("Insplit graph~~~~~~~~~")
  print_snf(gi)
  print("Homology")
  print_homology(g,gi)


def premade_graphs(og, ig):
  g = TwoGraph(og)
  gi = TwoGraph(ig)
  H0,H1,H2 = check_homology(g)
  print_generators("OG",H0,H1,H2)
  H0,H1,H2 = check_homology(gi)
  print_generators("\nIG",H0,H1,H2)


if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: <file name of original graph> <file name of insplit graph>")
    sys.exit(1)

  g = TwoGraph(sys.argv[1])
  gi = TwoGraph(sys.argv[2])

  print_homology(g,gi,cohomology=False)
  print()
  print_homology(g,gi,cohomology=True)

  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1,base_ring=ZZ)
  print("d1")
  print(d1)
  print("d2")
  print(d2)
  print("C.differential", C.differential(2))
  print("C.differential.dual", C.dual().differential(1))
