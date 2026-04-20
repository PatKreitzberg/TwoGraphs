from sympy import Matrix, pprint
from sympy.matrices.normalforms import smith_normal_form

def get_snf_generators(A):
    # A is your boundary matrix (m x n)
    # We use SymPy's smith_normal_form but need U and V
    # A more robust way to get U and V in SymPy is via the Smith Normal Form algorithm
    # or using the manual row/column operations.

    # For smaller matrices, we can compute the Smith Normal Form and the transforms:
    # Note: SymPy's smith_normal_form doesn't return U, V directly.
    # We can use the 'smith_normal_form' and then find U and V.

    # Alternatively, for generators in integer homology:
    D, U, V = A.smith_normal_form(with_transform=True)

    rank = A.rank()

    # 1. Image Generators:
    # Take the columns of U_inv corresponding to non-zero diagonal in D
    U_inv = U.inv()
    # The image is spanned by the first 'rank' columns of U_inv multiplied by diagonal elements
    # Since we want the basis for the image subgroup:
    img_generators = [U_inv.col(i) * D[i,i] for i in range(rank)]

    # 2. Kernel Generators:
    # Take the columns of V corresponding to the zero columns in D
    num_cols = A.cols
    ker_generators = [V.col(i) for i in range(rank, num_cols)]

    return D, img_generators, ker_generators


# Example matrix from your u, v graph
# Rows: u, v | Columns: r1, r2, b1, b2
A = Matrix([
    [0, 0,  1, -1],
    [0, 0, -1,  1]
])

pprint(A)

# 2. To get the spanning sets for Ker and Im specifically using SNF:
# We need U and V such that UAV = D
# In SymPy, we use the 'inv_normal_forms' or perform manual augmentation
# But for standard use, we can derive the spanning columns this way:

rank = A.rank()
num_columns = A.cols

# --- THE IMAGE ---
# The image is spanned by the first 'rank' columns of the matrix
# that would be U^-1. In practice, these are the columns of A
# that correspond to the pivot variables.
im_span = A.columnspace()

# --- THE KERNEL ---
# The kernel is spanned by the columns of V corresponding to the
# zero columns in D.
ker_span = A.nullspace()


print("\nColumns spanning the Image (im d1):")
for vec in im_span:
    print(vec)

print("\nColumns spanning the Kernel (ker d1):")
for vec in ker_span:
    print(vec)






# Just with SNF

D, img_gen, ker_gen = get_snf_generators(A)

print("Diagonal Matrix D:")
print(D)

print("\nGenerators of Image (im d1):")
for g in img_gen:
    print(g)

print("\nGenerators of Kernel (ker d1):")
for g in ker_gen:
    print(g)
