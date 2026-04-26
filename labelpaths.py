import numpy as np

def get_labeled_edges(matrix, label):
  """Returns a dictionary where keys are (i, j) and values are lists of labeled edges."""
  edges = {}
  n = len(matrix)
  for i in range(n):
    for j in range(n):
      count = int(matrix[i, j])
      if count > 0:
        # Label format: (type, start, end, k)
        edges[(i, j)] = [(label, i, j, k + 1) for k in range(count)]
  return edges

def solve_matrix_paths(R_mat, B_mat):
  n = len(R_mat)
  r_edges = get_labeled_edges(R_mat, 'R')
  b_edges = get_labeled_edges(B_mat, 'B')

  # Dictionary to store pairs: (i, j) -> [(path_RB, path_BR), ...]
  path_pairs = {}

  for i in range(n):
    for j in range(n):
      paths_RB = []
      paths_BR = []

      # Calculate RB paths: vi -> vk (R) -> vj (B)
      for k in range(n):
        edges_rk = r_edges.get((i, k), [])
        edges_kj = b_edges.get((k, j), [])
        for e1 in edges_rk:
          for e2 in edges_kj:
            paths_RB.append((e1, e2))

      # Calculate BR paths: vi -> vk (B) -> vj (R)
      for k in range(n):
        edges_bk = b_edges.get((i, k), [])
        edges_kj = r_edges.get((k, j), [])
        for e1 in edges_bk:
          for e2 in edges_kj:
            paths_BR.append((e1, e2))

      # Pair them off (assuming len(paths_RB) == len(paths_BR))
      path_pairs[(i, j)] = list(zip(paths_RB, paths_BR))

  return path_pairs


# --- Example Usage ---
R = np.array([[1, 1, 1], [0, 1, 1], [1, 0, 1]])
B = np.array([[2, 1, 1], [0, 2, 1], [1, 0, 2]])


R = np.array([[2, 1], [1, 1]])
B = np.array([[3,1],[1,2]])

R = np.array([[1, 0],[0,1]])
B = np.array([[0,2],[0,0]])


results = solve_matrix_paths(R, B)

for (start, end), pairs in results.items():
  if pairs:
    print(f"Paths from v{start} to v{end}:")
    for idx, (rb, br) in enumerate(pairs):
      print(f"  Pair {idx+1}:")
      print(f"  RB: {rb[0]} -> {rb[1]}")
      print(f"  BR: {br[0]} -> {br[1]}")

print("R")
print(R)
print("B")
print(B)
print("RB")
print(R*B)
print("BR")
print(B*R)
