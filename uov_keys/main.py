from functions import *

m = 6
n = 16
F = GF(4,"x")

O, O_I = generate_oil_subspace_matrix(F,m,n) # Vertical matrices

list_P = generate_public_matrices(F,m,n,O)
P_1 = list_P[0]
print(O_I.transpose()*P_1*O_I)
print("\n")
list_M = generate_list_M(list_P)
M_1 = list_M[0]
print(O_I.transpose()*M_1*O_I)

