def print_homology(H0,H1,H2, i_H0, i_H1, i_H2,is_cohomology):
  if is_cohomology:
    print(f"Cohomology Graph\tInsplit: \nH_0 {H0}\t\t{i_H0} \nH_1 {H1}\t\t{i_H1} \nH_2 {H2}\t\t{i_H2}")
  else:
    print(f"Homology Graph\tInsplit: \nH^0 {H0}\t\t{i_H0} \nH^1 {H1}\t\t{i_H1} \nH^2 {H2}\t\t{i_H2}")
  print()

def  print_adj_matrices(R, B):
  print()
  print("Red matrix")
  for r in R:
    print(r)
  print()

  print("Blue matrix")
  for r in B:
    print(r)
  print()
