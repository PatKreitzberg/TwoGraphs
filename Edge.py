class Edge:
  def __init__(self, label, s, r):
    # label is string name of edge
    # s is source vertex
    # r is range vertex

    self.label = label
    self.s = s # source of edge
    self.r = r # range of edge
    self.degree_index = None
    self.range_of_commuting_squares = set()
    self.commuting_squares = set()

  def __contains__(self, v):
    return (v == self.s) or (v == self.r)
  def __str__(self):
    return self.label
  def __eq__(self, other):
    return self.label == other.label
  def __hash__(self):
    return hash(self.label)
  def __lt__(self, other):
    return self.label < other.label

  def F(self, i,ell):
    assert i == 1
    assert ell in [0,1]
    if ell == 0:
      return self.r
    if ell == 1:
      return self.s
