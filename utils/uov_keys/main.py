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
n = 12
F = GF(4, 'x')

# Generate Oil subspace matrices Vertical O_I = [O] O_I^t * M * O_I = 0
#                                                [I]
O, O_I = generate_oil_subspace_matrix(F, m, n)

# Generate and check for full-rank matrices
M = regenerate_until_full_rank(F, m, n, O)

# Compute symplectic basis
from sage.matrix.symplectic_basis import symplectic_basis_over_field
L = []
for i in range(0,m):
    L.append(symplectic_basis_over_field(M[i])[1])

# Get the standard m-dimensional isotropic subspace basis of M[i] 
L_submatrix = [L[i][0:m] for i in range(m)]
print("Check L_submatrix[0]*M[0]*L_submatrix[0]^t = 0:", (L_submatrix[0]*M[0]*L_submatrix[0].transpose()).is_zero())

# Compute P such that P*M[0]*P^t = M[1]
Q = diagonalize_full_alternating_matrix(F,M[0]) # Q * M[0] * Q^t = anti-id
R = diagonalize_full_alternating_matrix(F,M[1]) # R * M[1] * R^t = anti-id
P = R.inverse()*Q
print("Check: P*M[0]*P^t == M[1] ", P*M[0]*P.transpose() == M[1])

# Compute T (m*m) such that T*L_0 = L_1*P with L_i of size m *n and P size n*n
# Compute the pseudoinverse of L0
L_0_pseudo_inv = L_submatrix[0].pseudoinverse()

# Compute T
T = (L_submatrix[1] * P) * L_0_pseudo_inv

# Display the result
print("T = \n",T, "\nCheck T*L[0] == L[1]*P :", T*L_submatrix[0] == L_submatrix[1]*P," rank T is full ? ",T.rank()==m)

# Pick random full rank matrix A and compute B such that T = A.inverse() * B
while True:
    A = random_matrix(F, m, m)  # You can also use RR or ZZ depending on your field
    if A.rank() == m:  # Check if the matrix is full rank
        break

# Compute B such that T = A.inverse() * B
A_inv = A.inverse()
B = A*T

# Compute L = B*L_0 and L = A*L_1*P
# Compute L using the first formula: L = B * L_0
L1 = B * L_submatrix[0]

# Compute L using the second formula: L = A * L_1 * P
L2 = A * L_submatrix[1] * P

# Display results
print("Check if it yields the same L ",L1 == L2)

# Check if x, y in span(L_1 rows) then x*M[0]*y = 0
# Get the row space of L1
row_space_L1 = L1.row_space()

# Perform the check for all pairs of x, y in the row space of L1
def check_condition(M0,M1, row_space):
    basis = row_space.basis()
    for x in basis:
        for y in basis:
            # Compute x * M0 * y^T
            result_1 = x * M0 * y.column()
            result_2 = x * M1 * y.column()
            if result_1 == 0 and result_2 == 0:
                return True  # If any result is not zero, return False
    return False  # If all results are zero, return True

# Check the condition
is_valid = check_condition(M[0], M[1],row_space_L1)

# Output the result
if is_valid:
    print("The condition x*M[0]*y = x*M[1]*y = 0 DOES hold for all x, y in the span of the rows of L1.")
else:
    print("The condition x*M[0]*y = x*M[1]*y = 0 DOES NOT hold for all x, y in the span of the rows of L1.")

print(L1, "\n", " of rank ", L1.rank())

# Check
print(O_I.transpose())
# Get the row space of O_I.transpose()
row_space_O_I = O_I.transpose().row_space()
is_valid = check_condition(M[0], M[1],row_space_O_I)
# Output the result
if is_valid:
    print("The condition x*M[0]*y = x*M[1]*y = 0 DOES hold for all x, y in the span of the rows of O_I.")
else:
    print("The condition x*M[0]*y = x*M[1]*y = 0 DOES NOT hold for all x, y in the span of the rows of O_I.")
