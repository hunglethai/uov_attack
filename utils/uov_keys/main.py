from utils.uov_keys.functions import *
from utils.totallyisotropicsubspace.functions import *

m = 4
n = 10
F = GF(16,"x")

O, O_I = generate_oil_subspace_matrix(F,m,n) # Vertical matrices

list_P = generate_public_matrices(F,m,n,O)
list_M = generate_list_M(list_P)
for i in range(m):
    print(list_M[i])
    print("\n")

