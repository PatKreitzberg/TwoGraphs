'''
Prompt:
I have a file that looks like:
<Start file>
R = [[1, 2, 0, 0, 0],[0, 1, 0, 0, 0],[0, 0, 1, 0, 0],[0, 0, 0, 1, 0],[2, 1, 0, 0, 1]]
B = [[2, 4, 0, 0, 0],[0, 2, 0, 0, 0],[0, 0, 2, 0, 0],[0, 0, 0, 2, 0],[4, 6, 0, 0, 2]],
insplit vertex = 0

Graph commuting squares:
(E(deg=2, (0, 1)#=0 E(deg=1, (0, 0)#=0 ~ E(deg=1, (0, 1)#=0 E(deg=2, (0, 0)#=0)
(E(deg=2, (0, 1)#=1 E(deg=1, (0, 0)#=0 ~ E(deg=1, (0, 1)#=1 E(deg=2, (0, 0)#=0)

Insplit graph commuting squares:
(E(deg=2, (0, 1)#=3^1 E(deg=1, (0, 0)#=0^2 ~ E(deg=1, (0, 1)#=1^1 E(deg=2, (0, 0)#=1^2)
(E(deg=2, (0, 1)#=0^2 E(deg=1, (4, 0)#=0 ~ E(deg=1, (0, 1)#=0^2 E(deg=2, (4, 0)#=0)
<End file>

Write a python function to take in a file and return
R: nxn matrix made using regulart lists
B: nxn matrix made using regulart lists
commuting_squares: list of CommutingSquare objects

The CommutingSquare class is written below where r1, s1, r2, s2 are from the Edge class. For each line after
"Graph commuting squares:" there is a commuting square written as <Edge 1> <Edge 2> ~ <Edge 3> <Edge 4>

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

class Edge:
  def __init__(self, label, s, r, degree=None):
    # label is string name of edge
    # s is source vertex
    # r is range vertex

    self.label = label
    self.s = s # source of edge
    self.r = r # range of edge
    self.degree = degree
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

'''

import re
import ast

from Edge import Edge
from CommutingSquare import CommutingSquare

def parse_graph_file(filepath):
  with open(filepath, 'r') as f:
    content = f.read()

  # 1. Parse Matrices R and B
  # Searches for "R = [...]" and "B = [...]"
  r_match = re.search(r'R\s*=\s*(\[\[.*?\]\])', content, re.DOTALL)
  b_match = re.search(r'B\s*=\s*(\[\[.*?\]\])', content, re.DOTALL)
  vertex_match = re.search(r'insplit vertex\s*=\s*(\d+)', content)[1]

  R = ast.literal_eval(r_match.group(1)) if r_match else []
  # Remove trailing comma if it exists in the raw string before literal_eval
  b_str = b_match.group(1) if b_match else "[]"
  B = ast.literal_eval(b_str)

  # 2. Parse Commuting Squares
  commuting_squares = []

  # We only care about the section after "Graph commuting squares:"
  # but before any other section (like "Insplit graph commuting squares:")
  section_match = re.search(
    r'Graph commuting squares:\s*(.*?)(?:\n\n|\n[A-Z]|$)',
    content,
    re.DOTALL
  )

  if section_match:
    square_lines = section_match.group(1).strip().split('\n')

    # Regex to find: E(deg=X, (S, R)#=ID
    # Matches: deg=2, (0, 1)#=0
    edge_regex = r'E\(deg=(\d+),\s*\((\d+),\s*(\d+)\)#=(\d+)'

    for line in square_lines:
      if not line.strip() or '~' not in line:
        continue

      # Find all 4 edge definitions in the line
      matches = re.findall(edge_regex, line)
      if len(matches) == 4:
        edges = []
        for m in matches:
          deg, s, r, eid = int(m[0]), int(m[1]), int(m[2]), m[3]
          # Reconstruct the label as it appears in the file for the Edge class
          label = f"E(deg={deg}, ({s}, {r})#={eid}"
          edges.append(Edge(label, s, r, degree=deg))

        # Assign edges based on: (r1 s1 ~ r2 s2)
        sq = CommutingSquare(edges[0], edges[1], edges[2], edges[3])
        commuting_squares.append(sq)

  return R, B, commuting_squares, vertex_match

def relabel_edge(e):
  e.label = re.sub(r'\s+', '', e.label)
  return e


R, B, commuting_squares, insplit_vertex = parse_graph_file("./result_diff_h1_with_torsion.txt")
edges = set()
for cs in commuting_squares:
  cs.s1 = relabel_edge(cs.s1)
  cs.r1 = relabel_edge(cs.r1)
  cs.s2 = relabel_edge(cs.s2)
  cs.r2 = relabel_edge(cs.r2)
  cs.label = '(' + cs.r1.label + ' ' + cs.s1.label + ' ~ ' + cs.r2.label + ' ' + cs.s2.label + ')'

  edges.add(cs.s1)
  edges.add(cs.r1)
  edges.add(cs.s2)
  edges.add(cs.r2)

print("# edges of the form <edge label> <source vertex> <range vertex>")
deg_1_edges = "1 "
deg_2_edges = "1 "
for edge in edges:
  print(edge.label, edge.s, edge.r)
  if edge.degree == 1:
    deg_1_edges += edge.label + " "
  elif edge.degree == 2:
    deg_2_edges += edge.label + " "

print()
print("# degrees  of the form <degree> <edge> <edge>...")
print(deg_1_edges)
print(deg_2_edges)

print()
print("# commuting squares of the form <edge label> <edge label> ~  <edge label> <edge label>")
for cs in commuting_squares:
  print(cs)


print()
print("#notes")
print("insplit-v ", insplit_vertex)
