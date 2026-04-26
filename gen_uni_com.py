
# Generate B using the result from the previous script
matrix_B = generate_commuting_matrix(result, Z)

# Verification
AB = np.matmul(result, matrix_B)
BA = np.matmul(matrix_B, result)
is_commute = np.array_equal(AB, BA)

print("Matrix B (Commutes with A):")
print(matrix_B)
print(f"\nDoes AB = BA? {is_commute}")
