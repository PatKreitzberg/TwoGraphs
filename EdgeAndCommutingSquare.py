class Edge:
  def __init__(self, e, s, r):
    self.edge_label = e
    self.s = s
    self.r = r
    self.degree_index = None
    self.degree_tuple = None
  def __str__(self):
    return str(self.edge_label) + ": " + str(self.s) + " -> " + str(self.r)
  def __eq__(self, other):
    return self.edge_label == other.edge_label
  def __hash__(self):
    return hash(self.edge_label)

  def F(self, i,ell):
    assert i == 1
    assert ell in [0,1]
    if ell == 0:
      return self.r
    if ell == 1:
      return self.s

class CommutingSquare:
  def __init__(self, ll, lr, rl, rr):
    self.lhs = (ll,lr)
    self.rhs = (rl,rr)

    assert self.lhs[0].degree_index == self.rhs[1].degree_index
    assert self.lhs[1].degree_index == self.rhs[0].degree_index

    self.degree_indices = set([ll.degree_index, lr.degree_index, rl.degree_index, rr.degree_index])
    self.degree_index_to_edges = {
      self.lhs[0].degree_index: [self.lhs[0], self.rhs[1]],
      self.rhs[0].degree_index: [self.lhs[1], self.rhs[0]],
    }

  def __str__(self):
    return str(self.lhs[0].edge_label) + " " +str(self.lhs[1].edge_label) + "~" + str(self.lhs[0].edge_label) + " " +str(self.lhs[1].edge_label)

  def __getitem__(self, i):
    assert i < 4
    if i <= 1:
      return self.lhs[i]
    return self.rhs[i-2]

  def F(self, i, ell):
    # i is degree index
    assert i in self.degree_indices
    assert ell in [0,1]
    return self.degree_index_to_edges[i][ell]
