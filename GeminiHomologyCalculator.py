from sympy import Matrix, pprint
from sympy.matrices.normalforms import smith_normal_form, hermite_normal_form

def gemini_calculate_homology(d_i, d_next):
    d_i = Matrix(d_i)
    d_next = Matrix(d_next)

    # 1. Basis for Ker(d_i)
    kernel_basis = d_i.nullspace()
    if not kernel_basis:
        return "Group is 0"

    Z = Matrix.hstack(*kernel_basis)

    # 2. Solve for M such that Z * M = d2
    # We use d2.cols to iterate through each boundary vector
    M_columns = []
    for i in range(d_next.cols):
        # Solve for the coefficients of the i-th column of d2 in terms of Z
        sol, params = Z.gauss_jordan_solve(d_next.col(i))
        # sol will contain the coefficients (the 'coordinates' in the kernel)
        M_columns.append(sol)

    M = Matrix.hstack(*M_columns)

    # 3. SNF of M
    D = smith_normal_form(M)

    # 4. Interpret Results
    diagonal = [D[i, i] for i in range(min(D.shape))]
    free_rank = len(kernel_basis) - len([x for x in diagonal if x != 0])
    torsion = [x for x in diagonal if x > 1]

    return f"H_i = Z^{free_rank} + " + " + ".join([f"Z/{t}Z" for t in torsion])




def gemini_integer_homology(d1, d2):
    # 1. Integer Kernel Basis via HNF

    m, n = d1.shape
    augmented = d1.T.row_join(Matrix.eye(n))
    hnf = hermite_normal_form(augmented)

    ker_basis_list = []
    for i in range(hnf.rows):
        if all(val == 0 for val in hnf[i, :m]):
            ker_basis_list.append(hnf[i, m:].T)

    if not ker_basis_list:
        return "0"

    Z = Matrix.hstack(*ker_basis_list)

    # 2. Rectangular Solve for M
    # We solve Z * M = d2 column by column
    M_cols = []
    for i in range(d2.cols):
        # solve for d2's i-th column in terms of the basis Z
        sol, params = Z.gauss_jordan_solve(d2.col(i))
        M_cols.append(sol)

    M = Matrix.hstack(*M_cols)

    # 3. SNF for the Quotient
    D = smith_normal_form(M)

    # 4. Interpret
    diag = [D[i, i] for i in range(min(D.shape)) if D[i, i] != 0]
    free_rank = Z.cols - len(diag)
    torsion = [int(t) for t in diag if t > 1]

    parts = []
    if free_rank > 0:
        parts.append(f"Z^{free_rank}" if free_rank > 1 else "Z")
    for t in torsion:
        parts.append(f"Z/{t}Z")

    return " \oplus ".join(parts) if parts else "0"
