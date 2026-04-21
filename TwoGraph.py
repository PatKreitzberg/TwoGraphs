import re, math
from collections import defaultdict as dd
from fractions import Fraction
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

from Edge import Edge
from CommutingSquare import CommutingSquare
from BoundaryFunctionMatrix import BoundaryFunctionMatrix

class TwoGraph:
  def __init__(self, path):
    # Parse graph
    vertices, edge_label_to_edge, commuting_squares = self.parse(path)

    # To help calculate matrices
    self.vertices = list(vertices)
    self.edges = list(edge_label_to_edge.values())
    self.commuting_squares = list(commuting_squares)

    self.vertex_to_index = {v:i for i,v in enumerate(self.vertices)}
    self.edge_to_index  = {e:i for i,e in enumerate(self.edges)}
    self.commuting_square_to_index = {cs:i for i,cs in enumerate(self.commuting_squares)}

    # BOUNDARY FUNCTION MATRICES
    self.d_1 = BoundaryFunctionMatrix(
      self,
      1, # r = 1 so going from edges to vertices
      self.edges,
      self.vertices,
      self.edge_to_index,
      self.vertex_to_index,
      print=True
    )
    print('\n\n')
    self.d_2 = BoundaryFunctionMatrix(
      self,
      2, # r = 1 so going from edges to vertices
      self.commuting_squares,
      self.edges,
      self.commuting_square_to_index,
      self.edge_to_index,
      print=True
    )

  def print_commuting_squares(self):
    # Print commuting squares
      for cs in self.commuting_squares:
        print("Commuting square:", cs)
        for i in [1,2]:
          for ell in [0,1]:
            print(f"F_{i}^{ell} = {cs.F(i,ell).label}" )

  def parse(self, file_path):
    edge_label_to_edge = {}
    commuting_squares = set()
    vertices = set()
    current_section = None
    with open(file_path, 'r') as f:
      for line in f:
        line = line.strip()

        # Skip empty lines
        if not line:
          continue

        # Detect section headers
        if line.startswith('#'):
          header = line.lower()
          if 'edges' in header:
            current_section = 'edges'
          elif 'degrees' in header:
            current_section = 'degrees'
          elif 'commuting squares' in header:
            current_section = 'commuting_squares'
          elif 'notes' in header:
            return vertices, edge_label_to_edge, commuting_squares
          continue # skip the line that contains '#'

        # Parse based on the current section
        if current_section == 'edges':
          vertices,edges = self.parse_edge(line, vertices, edge_label_to_edge)
        elif current_section == 'degrees':
          edge_label_to_edge = self.parse_degree(line, edge_label_to_edge)
        elif current_section == 'commuting_squares':
          # Format: label_a label_b = label_c label_d
          if '~' in line:
            commuting_squares = self.parse_commuting_square(line, commuting_squares, edge_label_to_edge)
    return vertices, edge_label_to_edge, commuting_squares

  def parse_edge(self, line, vertices, edge_label_to_edge):
    # Format: <label> <v1> <v2>
    parts = line.split()
    if len(parts) == 3:
      e,s,r = parts
      edge = Edge(e,s,r)
      vertices.add(s)
      vertices.add(r)
      assert not (e in edge_label_to_edge.keys())
      edge_label_to_edge[e] = edge
    return vertices, edge_label_to_edge

  def parse_degree(self, line, edge_label_to_edge):
    parts = line.split()
    if len(parts) > 1:
      degree_index = int(parts[0])
      for e in parts[1:]:
        edge_label_to_edge[e].degree_index = degree_index
    return edge_label_to_edge

  def parse_commuting_square(self, line, commuting_squares, edge_label_to_edge):
    left_side, right_side = line.split('~')
    left_edge_1,  left_edge_2    = left_side.strip().split()
    right_edge_1, right_edge_2 = right_side.strip().split()

    # ab ~ cd
    commuting_squares.add(CommutingSquare(
      edge_label_to_edge[left_edge_1],
      edge_label_to_edge[left_edge_2],
      edge_label_to_edge[right_edge_1],
      edge_label_to_edge[right_edge_2]))
    return commuting_squares
