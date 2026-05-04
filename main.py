import sys, random
import numpy as np
from argparse import ArgumentParser

from sage.all import matrix, ZZ, ChainComplex

from python.RandomlyGeneratedTwoGraphForInsplitting import RandomlyGeneratedTwoGraphForInsplitting
from python.InsplitTwoGraph import InsplitTwoGraph
from python.RandomlyGeneratedTwoGraph import RandomlyGeneratedTwoGraph
from python.TwoGraph import TwoGraph

def sage_cohomology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  cochain = ChainComplex({1: d1, 2: d2}, degree=-1,base_ring=ZZ).dual()
  C0 = cochain.homology(0)
  C1 = cochain.homology(1)
  C2 = cochain.homology(2)
  return C0,C1,C2

def calc_homology(g):
  d1 = matrix(ZZ, g.d_1.matrix)
  d2 = matrix(ZZ, g.d_2.matrix)
  C = ChainComplex({1: d1, 2: d2}, degree=-1)
  H0 = C.homology(0)
  H1 = C.homology(1)
  H2 = C.homology(2)
  return {'H0':H0, 'H1':H1, 'H2':H2}

def save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2, prepend=''):
    fname = './results/graphs/'
    fname += prepend
    fname += 'n_vertices=' + str(g.n)
    fname += ' GH1:' +str(g_H1)
    fname += ' GH2:' +str(g_H2)
    fname += ' GIH1:' +str(gi_H1)
    fname += ' GIH2' +str(gi_H2)

    fname += ' z=' + str(g.z)
    end = ' ' + str(np.random.randint(0,1000000)) + '.txt'
    og_name = fname +' Original graph ' +  end
    ig_name = fname + ' insplit graph ' + end

    if (g_H1 !=  gi_H1) and (g_H2 !=  gi_H2):
        print(f"Diff homologies! Saving to \n{og_name} \nand\n{ig_name}")
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

    elif g_H2 !=  gi_H2:
      print("Diff H2!")
    else:
      print("Same homology :(")



def gen_random_graph_and_calc_homology_and_insplit_homology(n, z, runs, symmetric):
  print("n=", n, "z=", z)
  for run in range(runs):
    print('run', run)
    g = RandomlyGeneratedTwoGraphForInsplitting(n,z, symmetric=symmetric)

    if not g.is_legit:
      print("No insplit vertex found!")
      continue

    gi = InsplitTwoGraph(g, g.v, g.E1, g.E2)

    g_homology = calc_homology(g)
    gi_homology = calc_homology(gi)
    g_H0,  g_H1,  g_H2 =  g_homology['H0'],  g_homology['H1'], g_homology['H2']
    gi_H0, gi_H1, gi_H2 = gi_homology['H0'], gi_homology['H1'], gi_homology['H2']


    if 'C' in str(gi_H1) or 'C' in str(g_H1):
      print(f"G homology: \nH0 {g_H0}, \nH1 {g_H1} \nH2 {g_H2} ")
      print(f"In homology:\nH0 {gi_H0}, \nH1 {gi_H1} \nH2 {gi_H2}")

      print("Red matrix")
      for r in g.R:
        print(r)
      print()
      print("Blue matrix")
      for r in g.B:
        print(r)
      print()

      if g_H1 != gi_H1:
        print("Has torsion and different homologies!")
        save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2)


      if 'C' in str(g_H1) and 'C' not in str(gi_H1):
        print("WOW torsion in g_H1 but not gi_H1!")
        save_to_file(g, g_H0, g_H1, g_H2, gi, gi_H0, gi_H1, gi_H2, prepend='torsion-in-og-not-insplit')
        exit()



def premade_graphs(self, og, ig):
  g = TwoGraph(og)
  gi = TwoGraph(ig)


if __name__ == "__main__":
  parser = ArgumentParser()
  subparsers = parser.add_subparsers(dest="command", required=True)

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

  args = parser.parse_args()
  if args.command == 'random':
    n = args.n
    z = args.z
    runs = args.runs
    symmetric = args.symmetric
    gen_random_graph_and_calc_homology_and_insplit_homology(n, z, runs, symmetric)

  elif args.command == 'files':
    og_file,insplit_file = args.files
    premade_graphs(og_file, insplit_file)
