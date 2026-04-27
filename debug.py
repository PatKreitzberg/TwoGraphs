import numpy as np

from TwoGraph import TwoGraph


def p(m):
  s = ''
  for r in m:
    for c in r:
      s += str(c) + ' '
    s += '\n'
  print(s)

fname = './graphs/EFGGGP_Figure3.txt'
g = TwoGraph(load_from=fname)


p(g.d_1.matrix)
print()
p(g.d_2.matrix)
print()


A = np.matrix(g.d_1.matrix)
B = np.matrix(g.d_2.matrix)
print()
print('AB=\n',np.matmul(A,B))



d1 = [
  [0,  0,  1,  0,  1, 0],
  [0,  0, -1,  0, -1, 0]
]
d2 = [
  [0,  1,  0,  0],
  [0,  0, -1,  0],
  [0, -1,  1,  0],
  [0, -1,  0,  0],
  [0,  1, -1,  0],
  [0,  0,  1,  0]
]

print(d1)
print('d1d2=\n',np.matmul(np.matrix(d1),np.matrix(d2)))
