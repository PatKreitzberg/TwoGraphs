"""
Generate two n×n integer matrices A and B such that AB = BA.

Strategy:
  - A is random with entries in [0, z].
  - B = c0*I + c1*A + c2*A^2 + ... (a polynomial in A with positive integer
    coefficients). Any polynomial of A commutes with A.
  - The constant term c0 is chosen large enough so that every entry of B is
    positive (adding c0*I to a matrix that commutes still commutes).
  - All arithmetic is done in exact Python integers to avoid float rounding.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Pure-Python exact integer matrix helpers
# ---------------------------------------------------------------------------

def _identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]

def _matmul(X, Y, n):
    return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]

def _matadd(X, Y, n):
    return [[X[i][j] + Y[i][j] for j in range(n)] for i in range(n)]

def _matscale(c, X, n):
    return [[c * X[i][j] for j in range(n)] for i in range(n)]

def _matmin(X):
    return min(X[i][j] for i in range(len(X)) for j in range(len(X[0])))


# ---------------------------------------------------------------------------
# Main generator
# ---------------------------------------------------------------------------

def generate_commuting_matrices(n, z, degree=2):
    """
    Generate n×n integer matrices A and B such that AB = BA.

    Args:
        n:      Matrix dimension.
        z:      Upper bound (inclusive) for entries of A (entries in [0, z]).
        degree: Degree of the polynomial used to construct B (default 2).

    Returns:
        A: n×n integer ndarray with entries in [0, z].
        B: n×n positive integer ndarray that commutes with A.
    """
    rng = np.random.default_rng()

    # --- Generate A ---
    A_np = rng.integers(0, z + 1, size=(n, n))
    A = A_np.tolist()

    # --- Build B = c1*A + c2*A^2 + ... (no constant term yet) ---
    # coeffs are random positive integers in [1, 10].
    coeffs = list(map(int, rng.integers(1, 11, size=degree)))

    B = [[0] * n for _ in range(n)]
    A_power = [row[:] for row in A]   # A^1
    for c in coeffs:
        B = _matadd(B, _matscale(c, A_power, n), n)
        A_power = _matmul(A_power, A, n)

    # --- Choose c0 so all entries of c0*I + B are positive ---
    # Adding c0*I commutes with everything, so commutativity is preserved.
    min_entry = _matmin(B)
    c0 = max(1, 1 - min_entry)

    I = _identity(n)
    B = _matadd(B, _matscale(c0, I, n), n)

    return A_np, np.array(B)


# ---------------------------------------------------------------------------
# Verification and display (exact integer arithmetic)
# ---------------------------------------------------------------------------

def verify_commutes(A, B):
    """Return True if AB == BA exactly (using Python big-int arithmetic)."""
    A_list = A.tolist()
    B_list = B.tolist()
    n = len(A_list)
    AB = _matmul(A_list, B_list, n)
    BA = _matmul(B_list, A_list, n)
    return AB == BA

def print_matrix(name, M):
    data = M.tolist()
    width = max(len(str(v)) for row in data for v in row)
    print(f"\n{name}:")
    for row in data:
        print("  [" + "  ".join(f"{v:{width}d}" for v in row) + "]")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate commuting integer matrices A and B.")
    parser.add_argument("-n", type=int, default=3,
                        help="Matrix dimension (default: 3)")
    parser.add_argument("-z", type=int, default=9,
                        help="Max entry value for A (default: 9)")
    parser.add_argument("--degree", type=int, default=2,
                        help="Polynomial degree for B (default: 2)")
    args = parser.parse_args()

    print(f"Generating {args.n}x{args.n} commuting matrices "
          f"(A entries in [0, {args.z}], poly degree {args.degree})...")

    A, B = generate_commuting_matrices(args.n, args.z, args.degree)
    n = args.n

    print_matrix("A", A)
    print_matrix("B", B)

    AB = np.array(_matmul(A.tolist(), B.tolist(), n))
    BA = np.array(_matmul(B.tolist(), A.tolist(), n))
    print_matrix("AB", AB)
    print_matrix("BA", BA)

    ok = verify_commutes(A, B)
    print(f"\nAB == BA: {ok}")
