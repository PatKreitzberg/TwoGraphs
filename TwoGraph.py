import re, math
from collections import defaultdict as dd
from fractions import Fraction
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

from Edge import Edge
from CommutingSquare import CommutingSquare
from BoundaryFunctionMatrix import BoundaryFunctionMatrix

class TwoGraph:
  def __init__(self, load_from=None):
    vertices, edge_label_to_edge, commuting_squares = None, None, None
    self.R_degree = 1
    self.B_degree = 2

    if type(load_from) is str: # parse from graph
      vertices, edge_label_to_edge, commuting_squares = self.parse(load_from)
      self.vertices = list(vertices)
      self.edges = list(edge_label_to_edge.values())
      self.commuting_squares = list(commuting_squares)
      self.calculate_boundary_matrices()
    else:
      pass

  def save_to_file(self, filename, insplit_v = None):
    '''Write to file which can be opened in this class.

    '''
    nl = '\n'
    with open(filename, 'w+') as f:
      f.write('# edges of the form <edge label> <source vertex> <range vertex>')
      deg_1_str = '1 '
      deg_2_str = '1 '
      for e in self.edges:
        f.write(e.label + ' ' + e.r + ' ' + e.s + nl)
        if e.degree == 1:
          deg_1_str += e.label + ' '
        if e.degree == 2:
          deg_2_str += e.label + ' '
      f.write('')
      f.write('# degrees  of the form <degree> <edge> <edge>...')
      f.write(deg_1_str + '\n' + deg_2_str + nl)
      f.write('# commuting squares of the form <edge label> <edge label> ~  <edge label> <edge label>')
      for cs in self.commuting_squares:
        f.write(cs.s1 +' '+ cs.r1 +' '+  cs.s2 +' '+  cs.r2 + nl)

      if insplit_v is not None:
        f.write("#notes" + nl + 'insplit-v  ' + str(insplit_v))
      f.close()

  def print_path_matrix(self, A, title):
    adj_str = '['
    for row in A:
      adj_str += '['
      for col in row:
        adj_str += str(col) + ', '
      adj_str = adj_str[:-2]
      adj_str += ']\n'
    adj_str = adj_str[:-1] + ']'
    print(title)
    print(adj_str)

  def print_adj_matrices(self, title):
    R = [[0]*self.n for i in range(self.n)]
    B = [[0]*self.n for i in range(self.n)]
    for edge in self.edges:
      if edge.degree == self.R_degree:
        R[edge.s][edge.r] += 1
      else:
        B[edge.s][edge.r] += 1

    self.print_path_matrix(R, "Red adjacency matrix for" + title)
    self.print_path_matrix(B, "Blue adjacency matrix for" + title)


  def range_inverse_of_vertex(self, v):
    assert len(self.edges) > 0
    assert v in self.vertices
    return {e  for e in self.edges if e.r == v}

  def source_inverse_of_vertex(self, v):
    assert len(self.edges) > 0
    assert v in self.vertices
    return {e  for e in self.edges if e.s == v}

  def calculate_boundary_matrices(self):
    # To help calculate matrices
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
      calc_ker=False,
      calc_img=False
    )
    self.d_2 = BoundaryFunctionMatrix(
      self,
      2, # r = 1 so going from edges to vertices
      self.commuting_squares,
      self.edges,
      self.commuting_square_to_index,
      self.edge_to_index,
      calc_ker=True,
      calc_img=False
    )

  def parse(self, file_path):
    edge_label_to_edge = {}
    commuting_squares = set()
    vertices = set()
    current_section = None
    with open(file_path, 'r') as f:
      for line in f:
        line = line.strip()

        # Skip empty lines
        if not line or '%' in line:
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
      degree = int(parts[0])
      for e in parts[1:]:
        edge_label_to_edge[e].degree = degree
    return edge_label_to_edge

  def parse_commuting_square(self, line, commuting_squares, edge_label_to_edge):
    left_side, right_side = line.split('~')
    left_range_edge_label,  left_source_edge_label    = left_side.strip().split()
    right_range_edge_label, right_source_edge_label = right_side.strip().split()

    left_range_edge = edge_label_to_edge[left_range_edge_label]
    left_source_edge = edge_label_to_edge[left_source_edge_label]
    right_range_edge = edge_label_to_edge[right_range_edge_label]
    right_source_edge = edge_label_to_edge[right_source_edge_label]

    # ab ~ cd
    cs = CommutingSquare(left_range_edge, left_source_edge, right_range_edge, right_source_edge)
    commuting_squares.add(cs)

    left_range_edge.range_of_commuting_squares.add(cs)
    right_range_edge.range_of_commuting_squares.add(cs)

    range_v = left_range_edge.r
    source_v = left_source_edge.s

    left_range_edge.commuting_squares.add(cs)
    left_source_edge.commuting_squares.add(cs)
    right_range_edge.commuting_squares.add(cs)
    right_source_edge.commuting_squares.add(cs)
    return commuting_squares

  def draw_graph(self):
    g = nx.MultiDiGraph()
    deg_index_to_color = [None, "red", "blue"] # since degree > 0
    edge_list = [(e.s,e.r, {'color':deg_index_to_color[e.degree], 'label':e.label}) for e in self.edges]
    g.add_edges_from(edge_list)
    self.draw_multidigraph(g)

  def draw_multidigraph(self, g):
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("#1a1a2e")
    ax.set_facecolor("#1a1a2e")

    # --- Layout ---
    pos = nx.spring_layout(g, seed=42)

    # Spread nodes a bit more if only a few
    if len(pos) <= 3:
        keys = list(pos.keys())
        if len(keys) == 1:
            pos[keys[0]] = np.array([0.0, 0.0])
        elif len(keys) == 2:
            pos[keys[0]] = np.array([-0.5, 0.0])
            pos[keys[1]] = np.array([0.5, 0.0])

    # --- Draw edges ---
    # Group edges by (u, v) pair so we can fan them out
    edge_groups = {}
    for u, v, data in g.edges(data=True):
        key = (u, v)
        edge_groups.setdefault(key, []).append(data)

    label_positions = []  # [(x, y, label, color)]

    for (u, v), edges in edge_groups.items():
        n = len(edges)
        is_loop = u == v

        for i, data in enumerate(edges):
            color = data.get("color", "white")
            label = data.get("label", "")

            if is_loop:
                # Draw a self-loop as a circular arc above the node
                x, y = pos[u]
                angle_offset = (i - (n - 1) / 2) * 0.35  # fan loops sideways
                loop_radius = 0.12 + i * 0.04

                theta = np.linspace(0 + angle_offset, 2 * np.pi + angle_offset, 200)
                lx = x + loop_radius * np.cos(theta)
                ly = y + loop_radius * 1.5 * np.sin(theta) + loop_radius * 1.2

                ax.plot(lx, ly, color=color, linewidth=1.8, zorder=2)

                # Arrow at the end of the loop
                dx = lx[-1] - lx[-2]
                dy = ly[-1] - ly[-2]
                ax.annotate(
                    "",
                    xy=(lx[-1], ly[-1]),
                    xytext=(lx[-1] - dx * 3, ly[-1] - dy * 3),
                    arrowprops=dict(
                        arrowstyle="-|>",
                        color=color,
                        lw=1.5,
                        mutation_scale=14,
                    ),
                    zorder=3,
                )

                # Label at the top of the loop
                mid_idx = len(lx) // 2
                label_positions.append((lx[mid_idx], ly[mid_idx] + 0.03, label, color))

            else:
                # Fan out multiple edges between same pair using arc curvature
                # Spread: center edge is straight-ish, others curve left/right
                spread = 0.25
                if n == 1:
                    rad = 0.15  # slight curve even for single edge
                else:
                    rad = -spread + i * (2 * spread / (n - 1)) if n > 1 else 0.0

                src = pos[u]
                dst = pos[v]

                # Draw curved arrow
                arrow = FancyArrowPatch(
                    posA=src,
                    posB=dst,
                    connectionstyle=f"arc3,rad={rad:.3f}",
                    arrowstyle="-|>",
                    color=color,
                    linewidth=1.8,
                    mutation_scale=16,
                    zorder=2,
                    shrinkA=12,
                    shrinkB=12,
                )
                ax.add_patch(arrow)

                # Label at midpoint of the arc
                mx = (src[0] + dst[0]) / 2
                my = (src[1] + dst[1]) / 2
                # Offset perpendicular to the edge
                dx = dst[0] - src[0]
                dy = dst[1] - src[1]
                length = np.sqrt(dx**2 + dy**2) or 1
                perp = np.array([-dy, dx]) / length
                offset = perp * rad * 0.7
                lx = mx + offset[0]
                ly = my + offset[1]
                label_positions.append((lx, ly, label, color))

    # --- Draw nodes ---
    node_radius = 0.01
    for node, (x, y) in pos.items():
        circle = plt.Circle(
            (x, y),
            node_radius,
            color="#e0e0ff",
            zorder=4,
            linewidth=2,
            ec="#aaaacc",
        )
        ax.add_patch(circle)
        ax.text(
            x,
            y,
            str(node),
            ha="center",
            va="center",
            fontsize=14,
            fontweight="bold",
            color="#1a1a2e",
            zorder=5,
            fontfamily="monospace",
        )

    # --- Draw edge labels ---
    for lx, ly, label, color in label_positions:
        ax.text(
            lx,
            ly,
            label,
            ha="center",
            va="center",
            fontsize=9,
            color=color,
            fontfamily="monospace",
            bbox=dict(
                boxstyle="round,pad=0.2",
                fc="#1a1a2e",
                ec=color,
                alpha=0.85,
                linewidth=1,
            ),
            zorder=6,
        )

    # --- Legend for edge colors ---
    seen_colors = {}
    for _, _, data in g.edges(data=True):
        c = data.get("color", "white")
        if c not in seen_colors:
            seen_colors[c] = c
    legend_patches = [
        mpatches.Patch(color=c, label=c) for c in seen_colors
    ]
    ax.legend(
        handles=legend_patches,
        loc="lower right",
        framealpha=0.3,
        facecolor="#1a1a2e",
        edgecolor="#555577",
        labelcolor="white",
        fontsize=10,
    )

    ax.set_title(
        "K-Graph",
        color="#ccccee",
        fontsize=14,
        pad=12,
        fontfamily="monospace",
    )

    plt.tight_layout()
    plt.savefig("multidigraph.png", dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    print("Saved to multidigraph.png")
    plt.show()
