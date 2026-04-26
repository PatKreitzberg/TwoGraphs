import numpy as np

def generate_commuting_matrices(n, z):
    """
    Generates two n x n integer matrices A and B such that AB = BA.
    Entries are roughly bounded by z.
    """

    # 1. Create random diagonal matrices (these always commute)
    d1_entries = np.random.randint(-z, z + 1, size=n)
    d2_entries = np.random.randint(-z, z + 1, size=n)

    D1 = np.diag(d1_entries)
    D2 = np.diag(d2_entries)

    # 2. Generate a random unimodular matrix P (det(P) = 1)
    # This ensures P_inv is also an integer matrix.
    P = np.eye(n, dtype=int)

    # Apply a few random shear operations to "mix" the identity matrix
    for _ in range(n * 2):
        i, j = np.random.choice(n, size=2, replace=False)
        factor = np.random.randint(-2, 3) # Keep factors small to avoid huge numbers
        P[i] += factor * P[j]

    # 3. Calculate P_inverse
    # Since det(P) = 1, the inverse will be all integers.
    P_inv = np.linalg.inv(P).astype(int)

    # 4. Transform the diagonal matrices
    # A = P * D1 * P_inv, B = P * D2 * P_inv
    A = P @ D1 @ P_inv
    B = P @ D2 @ P_inv

    return A, B

# --- Example Usage ---
n = 10
z = 2
A, B = generate_commuting_matrices(n, z)

print(f"Matrix A:\n{A}\n")
print(f"Matrix B:\n{B}\n")

# Verification
AB = A @ B
BA = B @ A

print("Is AB == BA?")
print(np.array_equal(AB, BA))

if np.array_equal(AB, BA):
    print("\nSuccess! The matrices commute.")
