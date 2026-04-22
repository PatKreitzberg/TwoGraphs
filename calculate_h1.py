import numpy as np
from hsnf import smith_normal_form

def calculate_h1_gemini(d1, d2):
    """
    Calculates H1 = ker(d1) / img(d2) using the hsnf package.
    d1 and d2 should be numpy arrays.
    """
    # 1. Compute SNF of d1 to find the basis of the kernel
    # D1 = L1 @ d1 @ R1
    D1, L1, R1 = smith_normal_form(d1)

    # The rank of d1 is the number of non-zero entries in D1
    rank_d1 = np.count_nonzero(np.diagonal(D1))

    # The basis for the kernel of d1 consists of the last (columns - rank)
    # columns of the R1 matrix.
    K = R1[:, rank_d1:]

    # 2. Express d2 in the basis of the kernel.
    # Since img(d2) is a subset of ker(d1), we solve K @ X = d2.
    # Because K is part of a unimodular matrix (R1), we can solve this
    # by taking the inverse of R1, or more simply:
    # R1_inv @ K = [0, I]^T, so R1_inv @ d2 will give us X in the bottom rows.
    R1_inv = np.round(np.linalg.inv(R1)).astype(int)
    X = R1_inv[rank_d1:, :] @ d2

    # 3. Compute the SNF of X to find the relations
    # D_h, L_h, R_h = smith_normal_form(X)
    D_h, _, _ = smith_normal_form(X)

    # 4. Extract group structure
    diagonal_elements = np.diagonal(D_h)

    # Torsion: diagonal elements > 1
    torsion = [int(d) for d in diagonal_elements if d > 1]

    # Free rank: size of kernel minus rank of d2
    # The rank of d2 is the number of non-zero elements in D_h
    rank_d2 = np.count_nonzero(diagonal_elements)
    free_rank = K.shape[1] - rank_d2

    # Formatting the output
    group_parts = []
    if free_rank > 0:
        group_parts.append(f"Z^{free_rank}" if free_rank > 1 else "Z")
    for t in torsion:
        group_parts.append(f"Z/{t}")

    result = " + ".join(group_parts) if group_parts else "0"
    return {"result":result, "free_rank": free_rank, "torsion": torsion}


def calculate_h1_claude(d1, d2):
    """
    Compute H1 = ker(d1) / im(d2) for a 2-graph chain complex.

    Parameters
    ----------
    d1 : numpy.ndarray  (shape: n0 x n1)  boundary map C1 -> C0
    d2 : numpy.ndarray  (shape: n1 x n2)  boundary map C2 -> C1

    Returns
    -------
    str : description of H1 as an abelian group
    """
    d1 = np.array(d1, dtype=int)
    d2 = np.array(d2, dtype=int)

    n1 = d1.shape[1]  # dimension of C1

    # --- Step 1: kernel of d1 ---
    D1, L1, R1 = smith_normal_form(d1)
    diag_d1 = [D1[i, i] for i in range(min(D1.shape))]
    rank_d1 = sum(1 for v in diag_d1 if v != 0)
    ker_dim = n1 - rank_d1  # dimension of ker(d1) as a free Z-module

    # --- Step 2: image of d2, restricted to ker(d1) ---
    # The SNF diagonal of d2 directly gives the invariant factors of im(d2)
    # within C1; those that land in the kernel contribute to H1's torsion.
    D2, L2, R2 = smith_normal_form(d2)
    diag_d2 = [D2[i, i] for i in range(min(D2.shape))]
    # Nonzero entries of D2 are the invariant factors of im(d2)
    invariant_factors = [v for v in diag_d2 if v != 0]
    rank_d2 = len(invariant_factors)

    # --- Step 3: assemble H1 ---
    # Free rank = dim(ker d1) - rank(d2)  [rank-nullity in the quotient]
    free_rank = ker_dim - rank_d2

    # Torsion summands: invariant factors > 1
    torsion = [abs(v) for v in invariant_factors if abs(v) > 1]

    # --- Step 4: format result ---
    if free_rank == 0 and not torsion:
        return "H1 = 0"

    parts = []
    if free_rank == 1:
        parts.append("Z")
    elif free_rank > 1:
        parts.append(f"Z^{free_rank}")

    for t in torsion:
        parts.append(f"Z/{t}Z")

    return "H1 = " + " ⊕ ".join(parts)
