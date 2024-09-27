from utils.uov_keys.functions import *
from utils.totallyisotropicsubspace.functions import *

def regenerate_until_full_rank(F, m, n, O):
    """
    Check if all matrices in list_M are full rank, otherwise regenerate list_M.
    
    Args:
    F : Finite Field (GF)
    m : Number of rows of the matrices
    n : Number of columns of the matrices
    O : Oil subspace matrix 
    
    Returns:
    list_M : List of full-rank matrices
    """
    while True:
        # Generate the public matrices
        list_P = generate_public_matrices(F, m, n, O)
        
        # Generate list_M from list_P
        list_M = generate_list_M(list_P)
        
        # Check if all matrices in list_M are full rank
        all_full_rank = True
        for M in list_M:
            if M.rank() < n:
                all_full_rank = False
                break
        
        # If all matrices are full rank, return list_M
        if all_full_rank:
            return list_M

# Parameters
m = 4
n = 8
F = GF(4, 'x')

# Generate Oil subspace matrices Vertical O_I = [O] O_I^t * M * O_I = 0
#                                                [I]
O, O_I = generate_oil_subspace_matrix(F, m, n)

# Generate and check for full-rank matrices
M = regenerate_until_full_rank(F, m, n, O)

# List of L^t*M*L == anti-identity        
L = [] 
for i in range(m):
    L_i = diagonalize_full_alternating_matrix(F,M[i]).transpose()
    L.append(L_i)

T = L[1]*L[0].inverse()
print(M[0] == T.transpose()*M[1]*T)
TO_I = T*O_I
print((O_I.transpose()*M[0]*O_I).is_zero())
print((TO_I.transpose()*M[1]*TO_I).is_zero())
print(is_invariant_subspace(F,O_I.transpose(),T))

