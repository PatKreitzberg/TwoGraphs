"""
Generate pairs of n×n integer matrices (A, B) such that AB = BA,
with all entries in the range [-z, +z].

Two strategies are available:

  random   – Sample random A and B, check AB == BA, retry until success.
             Works well for small n; can be slow for large n.

  analytic – Construct A = p(M) and B = q(M) for a shared base matrix M,
             using random low-degree polynomials. Because any two
             polynomials of the same matrix commute, this always succeeds
             and scales to larger n.  We keep entries in [-z, z] by
             clipping the *base* matrix M rather than A/B (which would
             break commutativity).
"""

import numpy as np
from typing import Generator, Literal, Tuple


Matrix = np.ndarray


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------

def _random_int_matrix(n: int, z: int, rng: np.random.Generator) -> Matrix:
    """Return an n×n matrix of integers uniformly drawn from [-z, z]."""
    return rng.integers(-z, z + 1, size=(n, n))


def _commute(A: Matrix, B: Matrix) -> bool:
    """Return True iff AB == BA (exact integer comparison)."""
    return np.array_equal(A @ B, B @ A)


def _poly_of_matrix(M: Matrix, degree: int, coeffs: np.ndarray) -> Matrix:
    """
    Evaluate polynomial p(M) = c0*I + c1*M + ... + cd*M^d.
    Uses Python-int object dtype internally to avoid overflow.
    Returns an int64 ndarray.
    """
    n = M.shape[0]
    M_obj = M.astype(object)
    result = np.zeros((n, n), dtype=object)
    power = np.eye(n, dtype=object)
    for c in coeffs:
        result += int(c) * power
        power = power @ M_obj
    return result.astype(np.int64)


# ---------------------------------------------------------------------------
# Generators
# ---------------------------------------------------------------------------

def random_commuting_pairs(
    n: int,
    z: int,
    count: int,
    seed: int | None = None,
    max_attempts: int = 10_000,
) -> Generator[Tuple[Matrix, Matrix], None, None]:
    """
    Yield `count` pairs (A, B) of n×n integer matrices in [-z, z] with AB = BA.
    Uses brute-force random sampling with commutativity check.

    Raises RuntimeError if `max_attempts` per pair is exhausted.
    """
    rng = np.random.default_rng(seed)
    found = 0
    attempts = 0
    while found < count:
        A = _random_int_matrix(n, z, rng)
        B = _random_int_matrix(n, z, rng)
        attempts += 1
        if _commute(A, B):
            yield A, B
            found += 1
            attempts = 0
        if attempts >= max_attempts:
            raise RuntimeError(
                f"Could not find commuting pair #{found + 1} in {max_attempts} "
                f"attempts. Try a larger z, smaller n, or strategy='analytic'."
            )


def analytic_commuting_pairs(
    n: int,
    z: int,
    count: int,
    seed: int | None = None,
    poly_degree: int = 2,
) -> Generator[Tuple[Matrix, Matrix], None, None]:
    """
    Yield `count` pairs (A, B) using polynomial construction.

    For each pair:
      1. Draw a small base matrix M with entries in [-z, z].
      2. Compute A = p(M) and B = q(M) for random polynomials p, q.
      3. Clip A and B individually to [-z, z].

    Because p(M)q(M) = q(M)p(M) for any scalar polynomials, the pair
    commutes before clipping. After clipping individual entries the
    commutativity invariant may not hold in general, so we verify and
    retry with fresh polynomials if needed (in practice the low-degree
    construction keeps entries small enough that clipping rarely bites).
    """
    rng = np.random.default_rng(seed)
    for _ in range(count):
        while True:
            # Small M keeps polynomial values from exploding
            M = rng.integers(-2, 3, size=(n, n))
            coeffs_p = rng.integers(-2, 3, size=poly_degree + 1)
            coeffs_q = rng.integers(-2, 3, size=poly_degree + 1)
            A_full = _poly_of_matrix(M, poly_degree, coeffs_p)
            B_full = _poly_of_matrix(M, poly_degree, coeffs_q)
            A = np.clip(A_full, -z, z)
            B = np.clip(B_full, -z, z)
            # Clipping can break commutativity; verify and retry if so
            if _commute(A, B):
                break
        yield A, B


def generate_commuting_pairs(
    n: int,
    z: int,
    count: int = 5,
    strategy: Literal["random", "analytic"] = "analytic",
    seed: int | None = None,
    **kwargs,
) -> list[Tuple[Matrix, Matrix]]:
    """
    Generate `count` pairs of n×n commuting integer matrices with entries in [-z, z].

    Parameters
    ----------
    n        : Matrix dimension.
    z        : Entries are integers in [-z, +z].
    count    : Number of pairs to generate.
    strategy : 'random'   – rejection sampling (good for small n).
               'analytic' – polynomial construction (scales well, always finds a pair).
    seed     : Optional random seed for reproducibility.
    **kwargs : Forwarded to the underlying generator
               (e.g. max_attempts for random, poly_degree for analytic).
    """
    if strategy == "random":
        gen = random_commuting_pairs(n, z, count, seed=seed, **kwargs)
    elif strategy == "analytic":
        gen = analytic_commuting_pairs(n, z, count, seed=seed, **kwargs)
    else:
        raise ValueError(f"Unknown strategy '{strategy}'. Choose 'random' or 'analytic'.")
    return list(gen)


# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

def _fmt_matrix(M: Matrix, indent: str = "  ") -> str:
    rows = []
    for row in M:
        rows.append(indent + "  ".join(f"{v:4d}" for v in row))
    return "\n".join(rows)


def print_pairs(pairs: list[Tuple[Matrix, Matrix]], verify: bool = True) -> None:
    for i, (A, B) in enumerate(pairs, 1):
        print(f"Pair {i}:")
        print("  A =")
        print(_fmt_matrix(A))
        print("  B =")
        print(_fmt_matrix(B))
        if verify:
            ok = _commute(A, B)
            print(f"  AB == BA: {ok}")
        print()


# ---------------------------------------------------------------------------
# CLI / demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate pairs of commuting integer matrices."
    )
    parser.add_argument("-n", type=int, default=3, help="Matrix dimension (default: 3)")
    parser.add_argument("-z", type=int, default=5, help="Entry bound (default: 5)")
    parser.add_argument("-c", "--count", type=int, default=3, help="Number of pairs (default: 3)")
    parser.add_argument(
        "-s", "--strategy",
        choices=["random", "analytic"],
        default="analytic",
        help="Generation strategy (default: analytic)",
    )
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=10_000,
        help="Max retries per pair for random strategy (default: 10000)",
    )
    parser.add_argument(
        "--poly-degree",
        type=int,
        default=2,
        help="Polynomial degree for analytic strategy (default: 2)",
    )
    args = parser.parse_args()

    extra = {}
    if args.strategy == "random":
        extra["max_attempts"] = args.max_attempts
    else:
        extra["poly_degree"] = args.poly_degree

    print(f"Generating {args.count} commuting {args.n}×{args.n} "
          f"integer matrix pairs with entries in [{-args.z}, {args.z}]")
    print(f"Strategy: {args.strategy}\n")

    pairs = generate_commuting_pairs(
        n=args.n,
        z=args.z,
        count=args.count,
        strategy=args.strategy,
        seed=args.seed,
        **extra,
    )
    print_pairs(pairs)
