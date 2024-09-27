from functions import *

F = GF(4,'x')
n = 6
# A = create_anti_identity_matrix(F,n//2)
A = random_matrix(F,n,n)
A = A + A.transpose()
print("\n A = \n", A, "\n of rank ",A.rank())

# # Set identity matrix of size n * n
# I = matrix.identity(F,n)

# # Brute force search isotropic matrices size m * n

# list_iso_matrices = brute_force_search_isotropic_matrices(A,F,n//2)

# number_iso_matrices = len(list_iso_matrices)
# print(number_iso_matrices)


L = diagonalize_full_alternating_matrix(F,A)
print("L = \n ")
print(L, " of rank ",L.rank())
if L != None:
    print("Check \n") 
    print(L*A*L.transpose())