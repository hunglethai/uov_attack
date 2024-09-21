from functions import *

F = GF(2,'x')
n = 8
A = create_anti_identity_matrix(F,n//2)
print(A)

# Set identity matrix of size n * n
I = matrix.identity(F,n)

# Brute force search isotropic matrices size m * n

list_iso_matrices = brute_force_search_isotropic_matrices(A,F,n//2)

number_iso_matrices = len(list_iso_matrices)
print(number_iso_matrices)