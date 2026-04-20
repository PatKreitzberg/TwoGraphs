"""
kernel_image.py
===============
Integer kernel and image bases for matrices with integer entries.

Both functions return integer bases using exact arithmetic — no floating
point, no fractions in the output.

Approach
--------
image_basis:
    Computes the Hermite Normal Form of A.  sympy's hermite_normal_form
    returns the HNF with zero columns dropped, so its columns form a
    ℤ-basis for im(A) ⊆ ℤ^m directly.

kernel_basis:
    Uses sympy's rational nullspace, then clears denominators by scaling
    each vector by the LCM of its entry denominators to get integer vectors.

Dependencies: sympy  (pip install sympy)
"""

import math
from fractions import Fraction
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form


def image_basis(A: list[list[int]]) -> list[list[int]]:
    """
    Return an integer basis for the image (column space) of matrix A.

    Parameters
    ----------
    A : 2D list of ints, shape (m, n).

    Returns
    -------
    List of column vectors (each a list of ints) forming a ℤ-basis for
    im(A) ⊆ ℤ^m.  Returns an empty list when A is the zero matrix.
    """
    M = Matrix(A)
    # HNF(M) is upper-triangular with zero columns dropped,
    # so its columns are exactly a ℤ-basis for im(A).
    H = hermite_normal_form(M)
    basis = []
    for j in range(H.cols):
        col = [int(H[i, j]) for i in range(H.rows)]
        if any(x != 0 for x in col):
            basis.append(col)
    return basis


def kernel_basis(A: list[list[int]]) -> list[list[int]]:
    """
    Return an integer basis for the kernel (null space) of matrix A.

    Parameters
    ----------
    A : 2D list of ints, shape (m, n).

    Returns
    -------
    List of column vectors (each a list of ints) forming a ℤ-basis for
    ker(A) ⊆ ℤ^n.  Returns an empty list when the kernel is trivial.
    """
    M = Matrix(A)
    rational_ns = M.nullspace()   # exact rational arithmetic

    basis = []
    for v in rational_ns:
        # Scale by LCM of denominators to clear fractions.
        fracs = [Fraction(str(entry)) for entry in v]
        denom_lcm = 1
        for f in fracs:
            denom_lcm = denom_lcm * f.denominator // math.gcd(denom_lcm, f.denominator)
        basis.append([int(f * denom_lcm) for f in fracs])
    return basis


# ---------------------------------------------------------------------------
# Demo / smoke tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    def check(label, A):
        M = Matrix(A)
        m, n = len(A), len(A[0])
        print(f"\n{'='*52}")
        print(f"  {label}  ({m}×{n})")
        print("  A =")
        for row in A:
            print("   ", row)

        img = image_basis(A)
        ker = kernel_basis(A)

        print(f"\n  Image basis (rank = {len(img)}):")
        for v in img:
            print("   ", v)

        print(f"\n  Kernel basis (nullity = {len(ker)}):")
        for v in ker:
            print("   ", v)

        rn_ok = len(img) + len(ker) == n
        print(f"\n  Rank–nullity: {len(img)} + {len(ker)} = {len(img)+len(ker)}  "
              f"[n={n}]  {'✓' if rn_ok else '✗ FAIL'}")

        for v in ker:
            r = list(M * Matrix(v))
            ok = all(x == 0 for x in r)
            print(f"  ker check  A @ {v} = {r}  {'✓' if ok else '✗ FAIL'}")

    # check("∂₁ of triangle", [[-1, 0, -1],
    #                          [ 1,-1,  0],
    #                            [ 0, 1,  1]])

    # check("2×3 example",    [[ 1, 2, 3],
    #                            [ 4, 5, 6]])

    # check("Full-rank 2×2",  [[ 1, 2],
    #                            [ 3, 5]])

    # check("Zero 2×3",       [[ 0, 0, 0],
    #                            [ 0, 0, 0]])
    check("My ex", [[0,0,-1,1],[0,0,1,-1]])
