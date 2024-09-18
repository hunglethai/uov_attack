from functions import *
from tqdm import tqdm

F = GF(2,'x')
n = 10
A = create_anti_identity_matrix(F,n//2)
print(A)

# Brute force search isotropic subspace
result = find_totally_isotropic_subspace(A, F, n//2)

if result:
    print("Found totally isotropic subspace with basis:")
    for vec in result:
        print(vec)
        
    # Generate a matrix from the found basis
    isotropic_matrix = generate_matrix_from_basis(F,result)
    print("Matrix generated from the basis of the isotropic subspace:")
    print(isotropic_matrix)
else:
    print("No totally isotropic subspace found.")


# Check
print("Check: \n")
print(isotropic_matrix.transpose()*A*isotropic_matrix)