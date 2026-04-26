import numpy as np

def generate_commuting_matrices(n, z):
    # 1. Generate matrix A with random integers between 0 and z
    A = np.random.randint(0, z + 1, size=(n, n))

    # 2. Generate matrix B such that AB = BA
    # We use the property that any matrix commutes with a polynomial of itself.
    # To ensure B has positive integer entries, we can use B = A + kI
    # where k is a positive integer, or simply B = A + np.eye(n, dtype=int)
    # If A already has zeros and we need B to be strictly positive,
    # we can use B = A^2 + A + I or a similar construction.

    # Simple construction: B = A + Identity
    # This ensures B has positive entries if A's diagonal is >= 0
    # To be safer for "strictly positive", we can add 1 to every element
    # and use the identity. However, B = A + I is the most standard approach.

    I = np.eye(n, dtype=int)
    B = A + I

    # If you need B to have NO zero entries at all:
    # B = np.dot(A, A) + A + np.ones((n,n), dtype=int)

    return A, B

# Parameters
n = 10  # Size
z = 10 # Max value

A, B = generate_commuting_matrices(n, z)

# Verification
AB = np.dot(A, B)
BA = np.dot(B, A)
is_commuting = np.array_equal(AB, BA)

print("Matrix A:\n", A)
print("\nMatrix B (A + I):\n", B)
print("\nAB == BA?", is_commuting)
