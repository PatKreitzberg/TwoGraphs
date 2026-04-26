import numpy as np

def generate_constrained_matrices(n, z):
    # 1. Generate A with det = 1 or -1
    # We use a lower triangular matrix with 1s on diagonal and small off-diagonals
    A = np.eye(n, dtype=int)

    # Fill lower triangle with small integers to keep growth controlled
    for i in range(n):
        for j in range(i):
            A[i, j] = np.random.randint(0, z // 2)

    # Apply a few row/column swaps or additions to "hide" the triangular structure
    # This keeps det(A) as 1 * 1 * ... * 1 = 1
    for _ in range(n):
        idx = np.random.permutation(n)
        A[idx[0]] += A[idx[1]]

    # 2. Generate B = A + kI
    # Find the minimum value in A to determine the smallest possible k
    min_a = np.min(A)
    # k must make all entries > 0. If min_a is 0, k=1. If min_a is -5, k=6.
    k = max(1, 1 - min_a)

    B = A + k * np.eye(n, dtype=int)

    return A, B

# Execute
n, z = 10, 10
A, B = generate_constrained_matrices(n, z)

print("Matrix A:")
print(A)
print(f"Det(A): {round(np.linalg.det(A.astype(float)))}")

print("\nMatrix B (A + kI):")
print(B)

# Verification
print("\nCommutativity check (AB - BA):")
print(np.dot(A, B) - np.dot(B, A))
print(sum(sum(np.dot(A, B) - np.dot(B, A))))
