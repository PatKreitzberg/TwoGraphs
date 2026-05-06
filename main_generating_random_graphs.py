import re
import numpy as np

from python.InsplitTwoGraph import InsplitTwoGraph
from python.RandomlyGeneratedTwoGraphForInsplitting import RandomlyGeneratedTwoGraphForInsplitting

from main_calculate_homology import homology, cohomology

def verbose_output(g, gi, g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2):
  print()
  print("Verbose output")
  try:
    print("Insplit at", gi.v)
  except:
    pass

  print(f"Partitions E1={[str(e) for e in g.E1]}\tE2={[str(e) for e in g.E2]}")

  # Print homology
  print_homology(g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, False)

  # Print cohomology
  print_homology(g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2, True)

  # Print the red and blue adjacency matrices
  print_adj_matrices(g.R, g.B)


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

def calc_homologies(g, gi, verbose, only_og_torsion=True, any_torsion=False, g_h1_rank_gt=False):
    # Get homology and cohomology
    g_H0,  g_H1,  g_H2 = homology(g)
    gi_H0, gi_H1, gi_H2 = homology(gi)
    g_C0,  g_C1,  g_C2 = cohomology(g)
    gi_C0, gi_C1, gi_C2 = cohomology(gi)

    # Get rank and torsion
    g_H1_rank,  g_H1_torsion = parse_homology_rank(g_H1)
    gi_H1_rank, gi_H1_torsion = parse_homology_rank(gi_H1)

    #if verbose:
    #  verbose_output(g, gi, g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2)

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


    if any_torsion:
      print("in any torsion")
      if g_H1_torsion > 0 or gi_H1_torsion > 0:
        if g_H1 != gi_H1:
          print()
          print("nonzero torsion and different H1")
          if verbose:
            verbose_output(g, gi, g_H0, g_H1, g_H2, gi_H0, gi_H1, gi_H2, g_C0, g_C1, g_C2, gi_C0, gi_C1, gi_C2)
          save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2, prepend='just_diff_h1')

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
    calc_homologies(g, gi, verbose, any_torsion=True)
