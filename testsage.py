import numpy as np
from sage.all import matrix,ZZ


# 1. Define your matrix as a list of lists
data = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]

# 2. Convert to a Sage Matrix
# We specify 'ZZ' for integers or 'QQ' for rationals
A = matrix(ZZ, data)

# 3. Compute the Image (Column Space)
# Note: In Sage, .image() returns the row space of the matrix.
# To get the standard column-based image, use the matrix as is
# or transpose if you think in terms of column vectors.
img = A.image()

# 4. Compute the Kernel (Right Null Space)
# This finds vectors x such that A*x = 0
ker = A.right_kernel()

print("Matrix:")
print(A)
print("\nBasis for the Image:")
print(img.basis())
print("\nBasis for the Kernel:")
print(ker.basis())
