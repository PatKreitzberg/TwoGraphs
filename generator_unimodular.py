import numpy as np
import random

def generate_commuting_matrix(A, z):
    """
    Generates a matrix B such that AB = BA.
    Uses the property that B = c0*I + c1*A commutes with A.
    """
    n = A.shape[0]
    identity = np.eye(n, dtype=int)

    # We try small coefficients to stay within the bounds [-z, z]
    # B = c0*I + c1*A
    # We iterate to find coefficients that satisfy the integer bounds
    possible_coeffs = [c for c in range(-z, z+1) if c != 0]
    random.shuffle(possible_coeffs)

    for c1 in possible_coeffs:
        for c0 in possible_coeffs:
            B = c0 * identity + c1 * A
            if np.all(np.abs(B) <= z):
                return B

    # Fallback: if no combination works, the Identity matrix
    # is the trivial commuting matrix.
    return identity

def generate_unimodular_matrix(n, z, steps=100):
    """
    Generates an n x n matrix with determinant +/- 1 and
    integer values within the range [-z, z].
    """
    # Start with the Identity matrix (Det = 1)
    # To allow for Det = -1, we can randomly flip the sign of one row
    matrix = np.eye(n, dtype=int)
    if random.choice([True, False]):
        matrix[0, 0] = -1

    count = 0
    while count < steps:
        # Create a copy to test the next state
        next_matrix = matrix.copy()

        # Choose two distinct rows: i and j
        i, j = random.sample(range(n), 2)

        # Choose a random multiplier k
        # Small multipliers help keep the values from exploding too fast
        k = random.choice([-1, 1])

        # Elementary Operation: Row_i = Row_i + k * Row_j
        # This operation preserves the determinant.
        next_matrix[i] = next_matrix[i] + k * next_matrix[j]

        # Check if all elements are still within the bounds [-z, z]
        if np.all(np.abs(next_matrix) <= z):
            matrix = next_matrix
            count += 1
        else:
            # If out of bounds, we skip this step and try a different operation
            # This is a basic rejection sampling approach
            continue

    return matrix

# Parameters
N = 10  # Dimension
Z = 3 # Bounds

result = generate_unimodular_matrix(N, Z)
det = int(round(np.linalg.det(result)))

print(f"Generated {N}x{N} Matrix:")
print(result)
print(f"\nDeterminant: {det}")


# Generate B using the result from the previous script
matrix_B = generate_commuting_matrix(result, Z)

# Verification
AB = np.matmul(result, matrix_B)
BA = np.matmul(matrix_B, result)
is_commute = np.array_equal(AB, BA)

print("Matrix B (Commutes with A):")
print(matrix_B)
print(f"\nDoes AB = BA? {is_commute}")
