from python.Edge import Edge

class CommutingSquare:
  def __init__(self, r1, s1, r2, s2):
    # commuting squares need to have same source and range
    assert r1.r == r2.r
    assert s1.s == s2.s

    self.__setattr__('alpha', r1)

    self.__setattr__('s1', s1)
    self.__setattr__('r1', r1)
    self.__setattr__('s2', s2)
    self.__setattr__('r2', r2)
    self.__setattr__('path1', (r1, s1))
    self.__setattr__('path2', (r2, s2))
    self.__setattr__('label', '' + r1.label + ' ' + s1.label + ' ~ ' + r2.label + ' ' + s2.label + '')
    self.__setattr__('latex_label', '\substack{' + r1.label + ' ' + s1.label + r' \\ \sim \\ ' + r2.label + ' ' + s2.label + '}')

    assert self.path1[0].degree == self.path2[1].degree
    assert self.path1[1].degree == self.path2[0].degree

    self.degree_indices = set([r1.degree, s1.degree, r2.degree, s2.degree])


  def __setattr__(self, name, value):
    if hasattr(self, name) and name!='alpha' and name !='latex_label':
      print("WARNING: Changed", name)
      raise AttributeError(f"{name} is immutable")
    super().__setattr__(name, value)

  def __lt__(self, other):
    if self.alpha.label < other.alpha.label:
      return True
    return False

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

  def latex_str(self):
    return latex_label

  def F(self, i, ell):
    # i is degree index
    assert i in self.degree_indices
    assert ell in [0,1]

    ell_index = (ell+1)%2 # Because F_i^0 should be the range edge

    if self.path1[ell_index].degree == i:
      return self.path1[ell_index]

    return self.path2[ell_index]
