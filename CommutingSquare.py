from Edge import Edge

class CommutingSquare:
  def __init__(self, ll, lr, rl, rr):
    # l and r meaning left and right
    self.lhs = (ll,lr)
    self.rhs = (rl,rr)
    self.label = '(' + ll.label + ' ' + lr.label + ' ~ ' + rl.label + ' ' + rr.label + ')'
    self.r = {ll, rl} # is a set
    self.s = {lr, rr} # is a set

    assert self.lhs[0].degree_index == self.rhs[1].degree_index
    assert self.lhs[1].degree_index == self.rhs[0].degree_index

    self.degree_indices = set([ll.degree_index, lr.degree_index, rl.degree_index, rr.degree_index])
    self.degree_index_to_edges = {
      self.lhs[0].degree_index: [self.lhs[0], self.rhs[1]],
      self.rhs[0].degree_index: [self.lhs[1], self.rhs[0]],
    }

  def __eq__(self, other):
    return (self.lhs == other.lhs) and (self.rhs == other.rhs)

  def __contains__(self, edge):
    assert type(edge) is Edge
    return (edge in self.lhs) or (edge in self.rhs)

  def __hash__(self):
    return hash(self.label)

  def __str__(self):
    return self.label

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
