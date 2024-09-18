from functions import *
from tqdm import tqdm

F = GF(2,'x')
n = 6
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

# Search for matrix L such that isotropic_matrix.transpose() * L.transpose() * A * L * isotropic_matrix == 0
list_L = find_L_for_condition(isotropic_matrix, A, F)
print(len(list_L))
