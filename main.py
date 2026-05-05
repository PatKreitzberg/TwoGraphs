import re, ast
import random
import numpy as np
from argparse import ArgumentParser

from sage.homology.chain_complex import ChainComplex
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ


from python.RandomlyGeneratedTwoGraphForInsplitting import RandomlyGeneratedTwoGraphForInsplitting
from python.InsplitTwoGraph import InsplitTwoGraph
from python.RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph
from python.TwoGraph import TwoGraph

import math

def tikzpicture(gfile):
    degree_to_color = {1: "red, dashed", 2: "blue"}
    nl = '\n'

    g = TwoGraph(gfile)
    out = nl + r'\begin{tikzpicture}[' + nl
    out += r'  vertex/.style = {circle, draw, minimum size=0.8cm},' + nl
    out += r'  every loop/.style = {stealth-},' + nl
    out += r'  thick, ->, >=stealth' + nl
    out += r']' + nl + nl

    # 1. Position nodes in a circle to avoid (0,0) overlap
    radius = 3  # Adjust size of the graph layout
    num_vertices = len(g.vertices)
    for i, v in enumerate(g.vertices):
        angle = i * (360 / num_vertices)
        out += fr'\node[vertex] (v{v}) at ({angle}:{radius}) {{$ {v} $}};' + nl

    out += nl

    # 2. Track edge counts to prevent overlapping paths
    # Key: tuple of sorted vertex IDs for edges, or single ID for loops
    edge_counts = {}

    for e in g.edges:
        # Create a unique key for the pair (or loop)
        pair = tuple(sorted((e.s, e.r)))
        edge_counts[pair] = edge_counts.get(pair, 0) + 1
        count = edge_counts[pair]

        out += fr'\path[{degree_to_color[e.degree]}] (v{e.s}) edge ['

        if e.r == e.s:
            # Shift loop angle for each subsequent loop on the same node
            # Start at 0 degrees and move by 45 degrees each time
            angle = (count - 1) * 45
            out += fr'loop, out={angle+30}, in={angle-30}, looseness=8'
        else:
            # Vary the bend: 15, -15, 45, -45, etc.
            # Flips side and increases intensity to "stack" multiple edges
            bend = 20 * ((count + 1) // 2) * (1 if count % 2 == 1 else -1)
            out += fr'bend left={bend}'

        out += fr'] node[black, font=\small, auto] {{{e.label}}} (v{e.r});' + nl

    out += r'\end{tikzpicture}'
    print(out)




def oldtikzpicture(gfile):
  degree_to_color = {1:"red, dashed", 2:"blue"}
  nl = '\n'

  g = TwoGraph(gfile)
  out = '\n'
  out += r'\begin{tikzpicture}[' + nl
  out +=  r'vertex/.style = {circle, draw, minimum size=0.8cm},' + nl
  out +=  r'every loop/.style = {stealth-}, % Adds arrows to the loops' + nl
  out +=  r'thick,' + nl
  out +=  r'->, >=stealth, % Arrow style' + nl
  out +=  r'node distance=4cm' + nl
  out +=  r'],' + nl

  for v in g.vertices:
    out += r'\node[vertex] (v' + str(v) + r') at (0,0) {$' + str(v)  + '$};'
    out += nl + nl

  for e in g.edges:
    out += r'\path[' + degree_to_color[e.degree] + r'] (v' + str(e.s) + ') edge ['
    if e.r == e.s:
      out += 'loop right]'
    else:
      out += 'bend left=15]'
    out += ' node[black, below] {' + str(e.label) + '} (v' + str(e.r) + ');'
    out += nl + nl
  out += r'\end{tikzpicture}'
  print(out)


def print_vec(index_to_item, elt, chain, i):
  s = ''
  for ind,val in enumerate(chain._vec[i]):
    if val != 0:
      edge = index_to_item[ind]
      if val >0 and len(s) > 0:
        s+=' + '
      if val == 1:
        s += str(edge)
      elif  val == -1:
        s += ' - ' + str(edge)
      else:
        s += str(val) + str(edge)
  print(f'{elt}: ' + s)


def print_vecs(index_to_item, i, C):
  reprs = C.homology(i, generators=True)
  print("Reprs", reprs)
  for elt, chain in reprs:
    s = ''
    for ind,val in enumerate(chain._vec[i]):
      if val != 0:
        edge = index_to_item[ind]
        if val >0 and len(s) > 0:
          s+=' + '
        if val == 1:
          s += str(edge)
        elif  val == -1:
          s += ' - ' + str(edge)
        else:
          s += str(val) + str(edge)
    print(f'{elt}: ' + s)


def print_homology_structures(d1, d2, vertex_to_index, edge_to_index, commuting_square_to_index):
    # Invert the dictionaries for name lookups
    idx_to_v = {i: v for v, i in vertex_to_index.items()}
    idx_to_e = {i: e for e, i in edge_to_index.items()}
    idx_to_s = {i: s for s, i in commuting_square_to_index.items()}

    def pretty_print_vector(v, lookup):
        """Converts vector indices to object names."""
        parts = []
        for i, val in enumerate(v):
            if val == 0: continue
            label = lookup.get(i, f"idx_{i}")
            prefix = f"{val}*" if val != 1 else ""
            if val == -1: prefix = "-"
            parts.append(f"{prefix}[{label}]")
        return " + ".join(parts).replace("+ -", "- ")

    # 1. Kernel of d1: The 1-Cycles (Z1)
    # These are combinations of edges that form closed loops.
    print("=== 1-CYCLES: Ker(d1) ===")
    z1 = d1.right_kernel()
    if z1.rank() == 0:
        print("None (No cycles found)")
    else:
        for i, basis_vec in enumerate(z1.basis()):
            print(f"Cycle {i}: {pretty_print_vector(basis_vec, idx_to_e)}")

    print("\n" + "="*30 + "\n")

    # 2. Image of d2: The 1-Boundaries (B1)
    # These are edge cycles that are "filled in" by commuting squares.
    print("=== 1-BOUNDARIES: Im(d2) ===")
    # Using column_space() because each column of d2 represents the
    # boundary of a commuting square in the edge space.
    b1 = d2.column_space()
    if b1.rank() == 0:
        print("None (No boundaries found)")
    else:
        for i, basis_vec in enumerate(b1.basis()):
            print(f"Boundary {i}: {pretty_print_vector(basis_vec, idx_to_e)}")

    print("\n" + "="*30 + "\n")

    # 3. Relationship (Torsion Check)
    # If a cycle exists in Ker(d1) but its multiple exists in Im(d2),
    # that is where your C2 comes from.
    print("=== TORSION IDENTIFICATION ===")
    print(f"Rank of Ker(d1) [Cycles]: {z1.rank()}")
    print(f"Rank of Im(d2)  [Boundaries]: {b1.rank()}")
    print(f"H1 Rank (Betti Number): {z1.rank() - b1.rank()}")


def inspect_homology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1)
  print_homology_structures(d1, d2, g.vertex_to_index, g.edge_to_index, g.commuting_square_to_index)


def cohomology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  cochain = ChainComplex({1: d1, 2: d2}, degree=-1,base_ring=ZZ).dual()
  C0 = cochain.homology(0)
  C1 = cochain.homology(1)
  C2 = cochain.homology(2)
  return C0,C1,C2

def homology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1, base_ring=ZZ)
  H0 = C.homology(0)
  H1 = C.homology(1)
  H2 = C.homology(2)
  return H0,H1,H2

def save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2, prepend=''):
    fname = './results/graphs/'
    fname += prepend
    fname += 'n_vertices=' + str(g.n)
    fname += ' GH1:' +str(g_H1)
    fname += ' GH2:' +str(g_H2)
    fname += ' GIH1:' +str(gi_H1)
    fname += ' GIH2' +str(gi_H2)

    fname += ' z=' + str(g.z)
    end = ' ' + str(np.random.randint(0,10000000)) + '.txt'
    og_name = fname +' Original graph ' +  end
    ig_name = fname + ' insplit graph ' + end

    og_name = og_name.replace(" ", "_")
    ig_name = ig_name.replace(" ", "_")


    print(f"Saving to \n{og_name} \nand\n{ig_name}")
    g.save_to_file(og_name)
    gi.save_to_file(ig_name)

    # Append notes
    with open(og_name,'a') as f:
      f.write('#notes\n')
      f.write(' GH1:' +str(g_H1) + "\n")
      f.write(' GH2:' +str(g_H2) + "\n")
      f.write(' GIH1:' +str(gi_H1) + "\n")
      f.write(' GIH2:' +str(gi_H2) + "\n")
      f.close()

def print_homology(H0,H1,H2, i_H0, i_H1, i_H2,is_cohomology):
  if is_cohomology:
    print(f"Cohomology Graph\tInsplit: \nH_0 {H0}\t\t{i_H0} \nH_1 {H1}\t\t{i_H1} \nH_2 {H2}\t\t{i_H2}")
  else:
    print(f"Homology Graph\tInsplit: \nH^0 {H0}\t\t{i_H0} \nH^1 {H1}\t\t{i_H1} \nH^2 {H2}\t\t{i_H2}")
  print()

def print_adj_matrices(R,B):
  print("Red matrix")
  for r in R:
    print(r)
  print()
  print("Blue matrix")
  for r in B:
    print(r)
  print()

def verbose_output(g, gi, g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2):
  print()
  print("Verbose output")
  try:
    print("Insplit at", gi.v)
  except:
    pass

  # Print homology
  print_homology(g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, False)

  # Print cohomology
  print_homology(g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2, True)

  # Print the red and blue adjacency matrices
  print_adj_matrices(g.R, g.B)



def parse_homology_rank(homology_obj):
  """
  Parses the string representation of a Sage homology group
  to extract the free rank and the torsion value.
  """
  # Convert the object to a string
  s = str(homology_obj)

  # --- 1. Calculate Free Rank ---
  free_rank = 0

  if "Z^" in s:
    # Matches "Z^12" -> 12
    match = re.search(r'Z\^(\d+)', s)
    if match:
      free_rank = int(match.group(1))
  elif "Z x Z x Z x Z" in s:
    free_rank = 4
  elif "Z x Z x Z" in s:
    free_rank = 3
  elif "Z x Z" in s:
    free_rank = 2
  elif "Z" in s:
    # Ensure we aren't just seeing a 'Z' inside 'Z^' or a 'C'
    # This covers the single 'Z' case
    free_rank = 1
  else:
    free_rank = 0

  # --- 2. Calculate Torsion ---
  # Matches "C2", "C3", etc.
  torsion_match = re.search(r'C(\d+)', s)
  torsion = int(torsion_match.group(1)) if torsion_match else 0
  return free_rank, torsion


def calc_homologies(g, gi, verbose, only_og_torsion=True, any_torsion=False, g_h1_rank_gt=False):
    # Get homology and cohomology
    g_H0,  g_H1,  g_H2 = homology(g)
    gi_H0, gi_H1, gi_H2 = homology(gi)
    g_C0,  g_C1,  g_C2 = cohomology(g)
    gi_C0, gi_C1, gi_C2 = cohomology(gi)

    # Get rank and torsion
    g_H1_rank,  g_H1_torsion = parse_homology_rank(g_H1)
    gi_H1_rank, gi_H1_torsion = parse_homology_rank(gi_H1)

    if verbose:
      verbose_output(g, gi, g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2)

    ## if any H1 has torsion!
    # if g_H1_torsion > 0 or gi_H1_torsion > 0:
    #   verbose_output(g, gi, g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2)
    #   n_zeros = 0
    #   for row in g.R:
    #     for v in row:
    #       if v==0:
    #         n_zeros +=1
    #   if n_zeros > 1:
    #     save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2, prepend='small_with_torsion')
    #     print_adj_matrices(g.R,g.B)
    #     print("Exiting because example of small with torsion")
    #
    # if g_h1_rank_gt:
    #   print("rank diff")
    #   if g_H1_rank > gi_H1_rank:
    #     # easy to find examples with n=3 or n=4, z=1
    #     print()
    #     print("graph H1 rank greater than insplit H1 rank!")
    #     verbose_output(g, gi, g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2)
    #     save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2, prepend='GH1_rank_gt_GIH1')
    #     print("Exiting because g_H1_rank > gi_H1_rank")
    #     exit()


    # if any_torsion:
    #   print("in any torsion")
    #   if g_H1_torsion > 0 or gi_H1_torsion > 0:
    #     if g_H1 != gi_H1:
    #       print()
    #       print("nonzero torsion and different H1")
    #       if verbose:
    #         verbose_output(g, gi, g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2)
    #       save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2, prepend='just_diff_h1')

    if only_og_torsion:
      if g_H1_torsion > 0 and gi_H1_torsion == 0:
        print()
        print("Diff; WOW torsion in g_H1 but not gi_H1!")
        print(f"g.E1={[str(e) for e in g.E1]} and g.E2={[str(e) for e in g.E2]}")
        verbose_output(g, gi, g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2)
        save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2, prepend='wow-torsion-in-og-not-insplit')
        print("exiting because WOW torsion in g_H1 but not gi_H1!")
        return


def gen_random_graph_and_calc_homology_and_insplit_homology(n, z, runs, symmetric, verbose, R=None, B=None, only_og_torsion=True):
  print("n=", n, "z=", z)
  for run in range(runs):
    if run%50 == 0:
      print('run', run)

    g = RandomlyGeneratedTwoGraphForInsplitting(n,z, symmetric=symmetric)

    if not g.is_legit:
      print_legit(g)
      continue

    gi = InsplitTwoGraph(g, g.v, g.E1, g.E2)
    calc_homologies(g, gi, verbose)


def premade_graphs(og, ig):
  g = TwoGraph(og)
  gi = TwoGraph(ig)
  g_H0,  g_H1,  g_H2 = homology(g)
  gi_H0, gi_H1, gi_H2 = homology(gi)
  g_C0,  g_C1,  g_C2 = cohomology(g)
  gi_C0, gi_C1, gi_C2 = cohomology(gi)

  print('~~~~~~ Homology of G and Insplit G ~~~~~~')
  # Print homology
  print_homology(g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, False)
  # Print cohomology
  print_homology(g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2, True)
  print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
  inspect_homology(g)
  inspect_homology(gi)


def parse_matrix(M):
  if M is None:
    return None
  matrix = ast.literal_eval(M)
  assert type(matrix) is list
  return matrix

def print_legit(g):
  try:
    if g.v is None:
      print("G has no insplit vertex")
      return
  except:
    print("G has no v")
  try:
    if len(g.commuting_squares)==0:
      print("G has no commuting squares")
      return
  except:
    print("G has no commuting squares")





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

    inspect_homology(g)
    calc_homologies(g, gi, True)
    print("returning early")
    exit()
