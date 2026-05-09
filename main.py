import re, ast
import random
import numpy as np
from argparse import ArgumentParser

from sage.all import vector, matrix, ZZ, ChainComplex, latex

from main_generating_random_graphs import gen_random_graph_and_calc_homology_and_insplit_homology, calc_homologies
from main_calculate_homology import homology, cohomology
from main_print import print_homology

from python.RandomlyGeneratedTwoGraphForInsplitting import RandomlyGeneratedTwoGraphForInsplitting
from python.InsplitTwoGraph import InsplitTwoGraph
from python.RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph
from python.TwoGraph import TwoGraph

import math


def tikzpicture(gfile):
    degree_to_color = {1: "red, dashed", 2: "blue           "}
    nl = '\n'
    g = TwoGraph(gfile)
    loop_angle_delta = 70
    out = nl + r'\begin{tikzpicture}[' + nl
    out += r'  vertex/.style = {circle, draw, minimum size=0.8cm},' + nl
    out += r'  thick, ->, >=stealth,' + nl
    out += r'  loop spread/.style={looseness=8, distance=2cm},' + nl

    # We store the angles inside a TikZ "directory" called /angles/
    for v in g.vertices:
        out += fr'  /angles/v{v}/.initial=0,' + nl

    out += r']' + nl + nl

    # 1. Position nodes
    radius = 3
    num_vertices = len(g.vertices)
    for i, v in enumerate(g.vertices):
        angle = i * (360 / num_vertices)
        out += fr'\node[vertex] (v{v}) at ({angle}:{radius}) {{$ {v} $}};' + nl

    out += nl

    # 2. Edges and Loops
    edge_counts = {}
    for e in g.edges:
        pair = tuple(sorted((e.s, e.r)))
        edge_counts[pair] = edge_counts.get(pair, 0) + 1
        count = edge_counts[pair]

        out += fr'\path[{degree_to_color[e.degree]}] (v{e.s}) edge ['

        if e.r == e.s:
            # PGF math needs \pgfkeysvalueof to read the variable
            offset = (count - 1) * loop_angle_delta
            out += (fr'out={{\pgfkeysvalueof{{/angles/v{e.s}}} + {offset + 30}}}, '
                    fr'in={{\pgfkeysvalueof{{/angles/v{e.s}}} + {offset - 30}}}, loop')
        else:
            bend = 20 * ((count + 1) // 2) * (1 if count % 2 == 1 else -1)
            out += fr'bend left={bend}'

        out += fr'] node[black, auto] {{${e.label}$}} (v{e.r});' + nl

    out += r'\end{tikzpicture}'
    print(out)


def print_homology_and_solve_torsion(C, d1, d2, v_map, e_map, s_map):
    """
    Analyzes H1 homology, prints basis for cycles and boundaries,
    and dynamically solves for the squares that bound torsion cycles.
    """
    from sage.modules.free_module_element import vector

    # Invert maps
    idx_to_v = {i: v for v, i in v_map.items()}
    idx_to_e = {i: e for e, i in e_map.items()}
    idx_to_s = {i: s for s, i in s_map.items()}

    def pretty_print(v, lookup):
        parts = []
        for i, val in enumerate(v):
            if val == 0: continue
            label = lookup.get(i, f"idx_{i}")
            symbol = f"{val}*" if abs(val) != 1 else ("-" if val == -1 else "")
            parts.append(f"{symbol}[{label}]")
        return " + ".join(parts).replace("+ -", "- ")

    # 1. 1-Cycles (Ker d1)
    print("=== 1-CYCLES (Kernel of d1) ===")
    for i, b in enumerate(d1.right_kernel().basis()):
        print(f"Cycle {i}: {pretty_print(b, idx_to_e)}")

    print("\n=== 1-BOUNDARIES (Image of d2 by Square) ===")
    for j, col in enumerate(d2.columns()):
        if not col.is_zero():
            print(f"Square [{idx_to_s[j]}] bounds: {pretty_print(col, idx_to_e)}")
        else:
            print(f"Square [{idx_to_s[j]}] bounds: 0")

    print("\n=== 1-COMMUTING SQUARES (Kernel of d2) ===")
    for i, b in enumerate(d2.right_kernel().basis()):
        print(f"Cycle {i}: {pretty_print(b, idx_to_s)}")

    # 3. Dynamic Torsion Analysis
    print("\n=== DYNAMIC TORSION ANALYSIS ===")
    # Get H1 homology with generators
    h1 = C.homology(1, generators=True)

    for group, chain in h1:
        # Check if the group is finite (torsion)
        if group.is_finite():
            order = group.order()
            gen_vec = chain._vec[1]
            print(f"\nFound Torsion Group: {group}")
            print(f"Generator Chain: {pretty_print(gen_vec, idx_to_e)}")

            # Solve d2 * x = order * gen_vec
            try:
                # We want to find which squares (x) produce 'order' copies of the generator
                x = d2.solve_right(order * gen_vec)
                print(f"The {order}x multiple of this cycle is bounded by these squares:")
                print(f"  {pretty_print(x, idx_to_s)}")
            except (ValueError, RuntimeError):
                print(f"  [Note] Could not find a unique square combination in ZZ.")

    # 4. Free Part Analysis
    free_rank = sum(1 for group, _ in h1 if not group.is_finite())
    print(f"\nSummary: H1 has Betti number {free_rank} and {len(h1)-free_rank} torsion components.")


def print_simplified_boundaries_with_sources(d2, e_map, s_map):
    """
    Finds a simplified basis for the image of d2 and identifies
    which combination of commuting squares generates each basis vector.
    """
    idx_to_e = {i: e for e, i in e_map.items()}
    idx_to_s = {i: s for s, i in s_map.items()}

    def pretty_print(v, lookup):
        parts = []
        for i, val in enumerate(v):
            if val == 0: continue
            label = lookup.get(i, f"idx_{i}")
            # Formatting: 1*x -> [x], -1*x -> -[x], 2*x -> 2*[x]
            symbol = f"{val}*" if abs(val) != 1 else ("-" if val == -1 else "")
            parts.append(f"{symbol}[{label}]")
        if not parts: return "0"
        return " + ".join(parts).replace("+ -", "- ")

    # Compute the Hermite Normal Form basis for the column space
    img_d2 = d2.column_space()
    simplified_basis = img_d2.basis()

    print(f"=== SIMPLIFIED BOUNDARIES AND THEIR SOURCES ===")
    if not simplified_basis:
        print("No boundaries found.")
        return

    for i, boundary_vec in enumerate(simplified_basis):
        # 1. Format the boundary string (the Edges)
        boundary_str = pretty_print(boundary_vec, idx_to_e)

        # 2. Solve d2 * x = boundary_vec to find the squares
        # x will be a vector where x[j] is the coefficient for Square_j
        try:
            source_vec = d2.solve_right(boundary_vec)
            source_str = pretty_print(source_vec, idx_to_s)
        except (ValueError, RuntimeError):
            source_str = "Unknown combination (Inconsistent)"

        print(f"B{i}: {boundary_str.ljust(30)} Source: {source_str}")

def analyze_h1_generators(d1, d2, index_to_edge, index_to_square):
    d1 = d1.change_ring(ZZ)
    d2 = d2.change_ring(ZZ)

    # 1. Compute Right Kernel (Cycles)
    K = d1.right_kernel()
    basis_ker = K.basis_matrix() # Original basis from Sage
    k = K.rank()

    # 2. Coordinate matrix A: basis_ker.T * A = d2
    A = basis_ker.transpose().solve_right(d2)

    # 3. Smith Normal Form: S = U * A * V
    S, U, V = A.smith_form()

    # 4. The "Smith Basis" for the Kernel
    # The rows of U * basis_ker provide the basis for ker(d1)
    # that corresponds directly to the rows of S.
    smith_basis_ker = U * basis_ker

    # 1. Check if the square boundaries are actually cycles
    boundaries = d2.columns()
    for i, b in enumerate(boundaries):
        if not (d1 * b).is_zero():
            print(f"Square index {i} ({index_to_square[i]}) is NOT a cycle!")
            print(f"d1 * boundary = {d1 * b}")

    # 2. Check if the Smith Basis cycles are actually in the kernel
    for i in range(smith_basis_ker.nrows()):
        cycle = smith_basis_ker.row(i)
        if not (d1 * cycle.column()).is_zero():
            print(f"Smith Generator {i} is NOT in the kernel!")


    print("### Analysis of Cycles and their Boundaries\n")

    # We iterate through the rank of the kernel
    for i in range(k):
        # The cycle itself
        current_cycle_vec = smith_basis_ker.row(i)
        involved_edges = [
            (index_to_edge[j], coeff)
            for j, coeff in enumerate(current_cycle_vec) if coeff != 0
        ]
        cycle_desc = " + ".join([f"({c})*[{ed}]" for ed, c in involved_edges])

        # Check if this cycle is canceled, torsion, or free
        # This is determined by the diagonal of S (only exists up to rank of A)
        if i < A.rank():
            invariant_factor = S[i, i]

            # Find which squares bound this cycle via V
            square_combination = V.column(i)
            involved_squares = [
                (index_to_square[j], coeff)
                for j, coeff in enumerate(square_combination) if coeff != 0
            ]
            square_desc = " + ".join([f"({c})*[{sq}]" for sq, c in involved_squares])

            if invariant_factor == 1:
                print(f"--- GENERATOR {i+1} (CANCELED) ---")
                print(f"Cycle:    {cycle_desc}")
                print(f"Killed by: {square_desc}\n")
            else:
                print(f"--- GENERATOR {i+1} (TORSION Z/{invariant_factor}Z) ---")
                print(f"Cycle:    {cycle_desc}")
                print(f"Boundary: {invariant_factor} * Cycle = {square_desc}\n")
        else:
            print(f"--- GENERATOR {i+1} (FREE / UNBOUNDED) ---")
            print(f"Cycle:    {cycle_desc}\n")




def calculate_and_print_h1(d1, d2):
    d1 = d1.change_ring(ZZ)
    d2 = d2.change_ring(ZZ)

    # Check boundary condition explicitly
    product = d1 * d2
    if not product.is_zero():
        print("Boundary condition failed! d1 * d2 = ")
        print(product)
        return

    # Use RIGHT kernel for standard homology (column vector convention)
    K = d1.right_kernel()
    basis_ker = K.basis_matrix() # Rows are the basis vectors
    k = K.rank()

    try:
        # Solve (basis_ker.T) * X = d2
        # This finds coordinates of d2's columns in terms of the kernel basis
        A = basis_ker.transpose().solve_right(d2)
    except ValueError:
        print("Columns of d2 are not in the right kernel of d1.")
        return

    S, U, V = A.smith_form()

    assert S == U*A*V

    # Extract Homology
    r = A.rank()
    free_rank = k - r
    diag = [S[i,i] for i in range(min(S.nrows(), S.ncols()))]
    torsion = [d for d in diag if d > 1]

    # LaTeX Formatting
    parts = []
    if free_rank > 0: parts.append(f"\\mathbb{{Z}}^{{{free_rank}}}")
    for t in torsion: parts.append(f"\\mathbb{{Z}}/{t}\\mathbb{{Z}}")
    h1_latex = " \\oplus ".join(parts) if parts else "0"

    print("### Smith Normal Form (S)")
    print(f"$$\nS = {latex(S)}\n$$")
    print("\n### Homology Group")
    print(f"$$H_1 \cong {h1_latex}$$")


def inspect_homology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1)
  print_homology_and_solve_torsion(C, d1, d2, g.vertex_to_index, g.edge_to_index, g.commuting_square_to_index)
  print_simplified_boundaries_with_sources(d2, g.edge_to_index, g.commuting_square_to_index)


def premade_graphs(og, ig):
  g = TwoGraph(og)
  g_H0,  g_H1,  g_H2 = homology(g)
  g_C0,  g_C1,  g_C2 = cohomology(g)

  gi = TwoGraph(ig)
  gi_H0, gi_H1, gi_H2 = homology(gi)
  gi_C0, gi_C1, gi_C2 = cohomology(gi)

  # print('~~~~~~ Homology of G and Insplit G ~~~~~~')
  # # Print homology
  # print_homology(g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, False)
  # # Print cohomology
  # print_homology(g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2, True)
  # print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
  # print('~~~~~ Inspecting the original graph ~~~~~~')
  # inspect_homology(g)
  # print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
  # print('~~~~~ Inspecting the insplit graph ~~~~~~~')
  # inspect_homology(gi)
  # print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

  print("Smith for insplit")
  gid1,gid2 = matrix(ZZ, gi.d_1.matrix), matrix(ZZ, gi.d_2.matrix)
  print("Product\n", gid1*gid2)
  print("Gid1\n", gid1)
  print("Gid2\n", gid2)
  # --- Preparation for your 2-graph data ---
  # index_to_edge = {0: 'e1', 1: 'e2', ...}
  # index_to_square = {0: 'sq_ab_cd', 1: 'sq_gh_ij', ...}
  analyze_h1_generators(gid1, gid2, gi.index_to_edge, gi.index_to_commuting_square)




def parse_matrix(M):
  if M is None:
    return None
  matrix = ast.literal_eval(M)
  assert type(matrix) is list
  return matrix


if __name__ == "__main__":
  parser = ArgumentParser()
  subparsers = parser.add_subparsers(dest="command", required=True)

  # Single graph, want to get tikz of graph
  tikz_parser = subparsers.add_parser('tikz', help='Parser for randomly generated graphs')
  tikz_parser.add_argument(dest='gfile', help="Generate random graphs")

  # parse for when given files
  file_parser = subparsers.add_parser('files', help='Parser for randomly generated graphs')
  file_parser.add_argument("files", nargs=2, help="Generate random graphs")

  # parse for randomly generated
  runs_help = 'Number of runs (generate matrices and calc homology)'
  s_help = 'Make the red and blue adjacency matrices the same'
  random_parser = subparsers.add_parser('random', help='Parser for randomly generated graphs')
  random_parser.add_argument("n",      type=int,                  help='Number verices')
  random_parser.add_argument("z",       type=int,                  help='Entries in adj. matrices will be in [0,z]')
  random_parser.add_argument("runs",  type=int, default=1, help=runs_help)
  random_parser.add_argument("-s", "--symmetric", dest='symmetric',  required=False, default=False, action='store_true', help=s_help)
  random_parser.add_argument("-v", "--verbose", dest='verbose',  required=False, default=False, action='store_true', help='be verbose')

  # What to look for
  random_parser.add_argument("--only-og-torsion", dest='only_og_torsion',  required=False, default=True, action='store_true', help='be verbose')
  random_parser.add_argument("--any-torsion-any-diff-homologies", dest='any_torsion_any_diff_homologies',  required=False, default=False, action='store_true', help='be verbose')
  random_parser.add_argument("--g_h1_rank_gt", dest='g_h1_rank_gt',  required=False, default=False, action='store_true', help='be verbose')

  matrices_parser = subparsers.add_parser('matrices', help='Parser for matricesly generated graphs')
  matrices_parser.add_argument("-R",  dest='R',  required=False, default=None, help='Red adjacency matrix')
  matrices_parser.add_argument("-B",  dest='B',  required=False, default=None, help='Blue adjacency matrix')
  matrices_parser.add_argument("-v", "--verbose", dest='verbose',  required=False, default=False, action='store_true', help='be verbose')
  matrices_parser.add_argument("-vtx", "--vertex", dest='v',  required=False, type=int, help='be verbose')


  args = parser.parse_args()
  if args.command == 'random':
    n = args.n
    z = args.z
    runs = args.runs
    symmetric = args.symmetric
    verbose = args.verbose
    only_og_torsion = args.only_og_torsion
    any_torsion = args.any_torsion_any_diff_homologies
    g_h1_rank_gt = args.g_h1_rank_gt
    gen_random_graph_and_calc_homology_and_insplit_homology(n, z, runs, symmetric, verbose, R=None, B=None, only_og_torsion=only_og_torsion)

  elif args.command == 'tikz':
    gfile = args.gfile
    tikzpicture(gfile)

  elif args.command == 'files':
    og_file,insplit_file = args.files
    premade_graphs(og_file, insplit_file)

  elif args.command == 'matrices':
    for run  in range(100):
      print("arg.R is", args.R)
      print("arg.B is", args.B)
      verbose = args.verbose
      R = parse_matrix(args.R)
      B = parse_matrix(args.B)
      assert (R is None and B is None) or (type(R) is list and type(B) is list)

      n = len(R)
      assert len(R) == len(B)
      for i in range(len(R)):
        assert len(R[i]) == len(B[i])
        assert len(R[i]) == n

      Rmat = matrix(R)
      Bmat = matrix(B)
      assert Rmat*Bmat == Bmat*Rmat

      z = 0
      for mat in (R,B):
        for row in mat:
          rm = max(row)
          z = max(rm,z)
      assert z > 0

      print(f"Creating graph with R={R} and B={B}")
      g = RandomlyGeneratedTwoGraphForInsplitting(n,z,R=R,B=B)
      if not g.is_legit:
        print_legit(g)
        exit()

      gi = InsplitTwoGraph(g, g.v, g.E1, g.E2)
      calc_homologies(g, gi, True)
      print("Inspecting OG")
      inspect_homology(g)
      print("Inspecting Insplit")
      inspect_homology(gi)

      print("G edges")
      for e in g.edges:
        print(f'{e}: {e.s} -> {e.r}')
      print("GI edges")
      for e in gi.edges:
        print(f'{e}: {e.s} -> {e.r}')
      g_H0, g_H1, g_H2 = homology(g)
      gi_H0, gi_H1, gi_H2 = homology(gi)
      print_homology(g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, False)
