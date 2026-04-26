from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form

from sympy import Matrix

def claude_calculate_homology(d1: list[list[int]], d2: list[list[int]]):
    D1 = Matrix(d1)
    D2 = Matrix(d2)

    assert D1 * D2 == Matrix.zeros(D1.rows, D2.cols), "∂₁ ∘ ∂₂ ≠ 0"

    # --- H2 = ker(∂₂) ---
    rank_d2 = D2.rank()
    free_rank_H2 = D2.cols - rank_d2

    # --- H1 = ker(∂₁) / im(∂₂) ---
    # Step 1: integer basis for ker(∂₁)
    ker_basis = D1.nullspace()          # vectors over ℚ; clear denominators
    K_cols = [v * lcm_denominators(v) for v in ker_basis]
    K = Matrix.hstack(*K_cols)          # columns span ker(∂₁) over ℤ

    # Step 2: express im(∂₂) in this basis by solving K @ A = D2
    # Since im(D2) ⊆ ker(D1) by the chain condition, this has exact integer solutions
    A, _ = K.solve(D2, method='iQR')   # or solve column by column

    # Step 3: SNF of A gives invariant factors of the quotient
    A_snf = A.smith_normal_form()
    diagonal = [A_snf[i, i] for i in range(min(A_snf.rows, A_snf.cols))]

    torsion_H1  = [d for d in diagonal if d > 1]
    free_rank_H1 = K.cols - rank_d2

    def group_str(free_rank, torsion):
        parts = []
        if free_rank > 0:
            parts.append("ℤ" if free_rank == 1 else f"ℤ^{free_rank}")
        parts += [f"ℤ/{d}ℤ" for d in torsion]
        return " ⊕ ".join(parts) if parts else "0"

    return {
        "H1": {"free_rank": free_rank_H1, "torsion": torsion_H1,
               "group": group_str(free_rank_H1, torsion_H1)},
        "H2": {"free_rank": free_rank_H2, "torsion": [],
               "group": group_str(free_rank_H2, [])}
    }

def lcm_denominators(v):
    from math import lcm
    from functools import reduce
    denoms = [x.q for x in v if hasattr(x, 'q')]
    return reduce(lcm, denoms, 1)

# def claude_calculate_homology(d1: list[list[int]], d2: list[list[int]]):
#     """
#     Compute H1 = ker(d1) / im(d2) and H2 = ker(d2).

#     Args:
#         d1: boundary map ∂₁ as a list of lists of ints (shape: C0 x C1)
#         d2: boundary map ∂₂ as a list of lists of ints (shape: C1 x C2)

#     Returns:
#         dict with keys 'H1' and 'H2', each a dict describing the group.
#     """
#     D1 = Matrix(d1)
#     D2 = Matrix(d2)


#     # Sanity check
#     assert D1 * D2 == Matrix.zeros(D1.rows, D2.cols), "∂₁ ∘ ∂₂ ≠ 0 — not a valid chain complex"

#     # --- H2 = ker(∂₂) ---
#     # SNF of D2; invariant factors tell us rank of image
#     # ker has dimension = cols(D2) - rank(D2)
#     _, pivot_cols = D2.rref()  # over ℚ, but rank is the same over ℤ
#     rank_d2 = len(pivot_cols)
#     free_rank_H2 = D2.cols - rank_d2
#     # H2 is always free (it's a subgroup of a free abelian group with no quotient)

#     # --- H1 = ker(∂₁) / im(∂₂) via SNF of D2 restricted to ker(∂₁) ---
#     # Equivalently: compute SNF of D2 directly — the invariant factors give H1.
#     # The quotient ker(∂₁)/im(∂₂) is computed by finding the SNF of D2
#     # viewed as a map into C1. The nonzero diagonal entries d_i give:
#     #   - d_i = 1  → trivial contribution
#     #   - d_i > 1  → torsion summand ℤ/d_iℤ
#     #   - zero rows in SNF beyond rank → free ℤ summands
#     #
#     # We use the Smith Normal Form of D2 directly.
#     D2_snf = smith_normal_form(D2)
#     diagonal = [D2_snf[i, i] for i in range(min(D2_snf.rows, D2_snf.cols))]

#     torsion_H1 = [d for d in diagonal if d > 1]
#     num_trivial = sum(1 for d in diagonal if d == 1)

#     # rank of ker(∂₁)
#     _, pivot_cols_d1 = D1.rref()
#     rank_d1 = len(pivot_cols_d1)
#     ker_d1_rank = D1.cols - rank_d1

#     # free rank of H1 = dim ker(∂₁) - rank im(∂₂)
#     free_rank_H1 = ker_d1_rank - rank_d2

#     def group_str(free_rank, torsion):
#         parts = []
#         if free_rank > 0:
#             parts.append(f"ℤ^{free_rank}" if free_rank > 1 else "ℤ")
#         parts += [f"ℤ/{d}ℤ" for d in torsion]
#         return " ⊕ ".join(parts) if parts else "0"

#     H1 = {"free_rank": free_rank_H1, "torsion": torsion_H1, "group": group_str(free_rank_H1, torsion_H1)}
#     H2 = {"free_rank": free_rank_H2, "torsion": [],           "group": group_str(free_rank_H2, [])}

#     return {"H1": H1, "H2": H2}
