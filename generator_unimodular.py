import numpy as np
import random

def generate_unimodular_matrix_positive(n, z, steps=200):
    """
    Generates an n x n matrix with det +/- 1 and
    integer values within the range [0, z].
    """
    # Start with Identity (Det = 1)
    matrix = np.eye(n, dtype=int)

    # Optional: Swap two rows to make Det = -1
    if random.choice([True, False]):
        matrix[[0, 1]] = matrix[[1, 0]]

    count = 0
    while count < steps:
        next_matrix = matrix.copy()

        # Pick two rows and a multiplier
        i, j = random.sample(range(n), 2)
        # Use a small k (-1 or 1) to navigate the space without hitting bounds too fast
        k = random.choice([-1, 1])

        # Elementary Operation: Ri = Ri + k * Rj
        next_matrix[i] = next_matrix[i] + k * next_matrix[j]

        # Check if all elements are within [0, z]
        if np.all(next_matrix >= 0) and np.all(next_matrix <= z):
            matrix = next_matrix
            count += 1
        else:
            # If we go out of bounds, try a different operation
            continue

    return matrix

def generate_commuting_matrix_positive(A, z):
    """
    Generates a matrix B such that AB = BA with entries in [0, z].
    """
    n = A.shape[0]
    identity = np.eye(n, dtype=int)

    # Try combinations of B = c0*I + c1*A
    # We look for non-negative coefficients to keep B entries positive
    for c1 in range(z + 1):
        for c0 in range(z + 1):
            if c0 == 0 and c1 == 0: continue

            B = c0 * identity + c1 * A
            if np.all(B >= 0) and np.all(B <= z):
                return B

    return identity # Fallback to Identity if no other simple combination fits

# Configuration
N = 10
Z = 20

# Generate A
matrix_A = generate_unimodular_matrix_positive(N, Z)
det_A = int(round(np.linalg.det(matrix_A)))

# Generate B
matrix_B = generate_commuting_matrix_positive(matrix_A, Z)

# Verification
AB = np.matmul(matrix_A, matrix_B)
BA = np.matmul(matrix_B, matrix_A)
is_commute = np.array_equal(AB, BA)

print(f"Matrix A (Det = {det_A}):")
print(matrix_A)
print(f"\nMatrix B (Commutes with A):")
print(matrix_B)
print(f"\nVerification: Does AB = BA? {is_commute}")
print(f"Check Bounds: Max val in A={matrix_A.max()}, Max val in B={matrix_B.max()}")
