from Edge import Edge

class CommutingSquare:
  def __init__(self, r1, s1, r2, s2):
    # commuting squares need to have same source and range
    assert r1.r == r2.r
    assert s1.s == s2.s

    self.path1 = (r1,s1)
    self.path2 = (r2,s2)
    self.s1 = s1
    self.r1 = r1
    self.s2 = s2
    self.r2 = r2
    self.label = '(' + r1.label + ' ' + s1.label + ' ~ ' + r2.label + ' ' + s2.label + ')'



    assert self.path1[0].degree == self.path2[1].degree
    assert self.path1[1].degree == self.path2[0].degree

    self.degree_indices = set([r1.degree, s1.degree, r2.degree, s2.degree])
    self.degree_to_edges = {
      self.path1[0].degree: [self.path1[0], self.path2[1]],
      self.path2[0].degree: [self.path1[1], self.path2[0]],
    }

  def __eq__(self, other):
    return (self.path1 == other.path1) and (self.path2 == other.path2)

  def __contains__(self, edge):
    assert type(edge) is Edge
    return (edge in self.path1) or (edge in self.path2)

  def __hash__(self):
    return hash(self.label)

  def __str__(self):
    return self.label

  def __getitem__(self, i):
    assert i < 4
    if i <= 1:
      return self.path1[i]
    return self.path2[i-2]

  def F(self, i, ell):
    # i is degree index
    assert i in self.degree_indices
    assert ell in [0,1]
    return self.degree_to_edges[i][ell]
