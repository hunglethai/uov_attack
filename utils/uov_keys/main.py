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
m = 4  # i.e o
n = 10 # i.e. v
F = GF(4, 'x')

# Generate Oil subspace matrices Vertical O_I = [O] O_I^t * M * O_I = 0
#                                                [I]
O, O_I = generate_oil_subspace_matrix(F, m, n)

# Generate and check for full-rank matrices
M = regenerate_until_full_rank(F, m, n, O)

# Compute symplectic basis
from sage.matrix.symplectic_basis import symplectic_basis_over_field
L = []
for i in tqdm(range(0,m), ncols = 100, desc = "Compute symplectic basis ... "):
    L.append(symplectic_basis_over_field(M[i])[1])

# Get a full-rank m-dimensional isotropic subspace basis of M[i] 
# Iterate over the starting row index to create submatrices of size m x m
L_submatrix = []
for j in range(m):
    submatrix_list = []
    for i in range(n - m + 1):
        # Extract submatrix with columns from i to i+m
        submatrix = L[j][i:i+m]
        # Check if the submatrix is full rank
        if submatrix.rank() == m:
            submatrix_list.append(submatrix)
    L_submatrix.append(submatrix_list)

# Check if they are indeed isotropic
all_check_pass = True
for i in range(m):
    for j in range(len(L_submatrix[i])):
        if L_submatrix[i][j]*M[i]*L_submatrix[i][j].transpose() != Matrix(F,m,m) or L_submatrix[i][j].rank() != m:
            all_check_pass = False
print("Check L_submatrix[i][j]*M[i]*L_submatrix[i][j]^t = 0: ", all_check_pass)

# Compute P such that P*M[0]*P^t = M[1]
Q = diagonalize_full_alternating_matrix(F,M[0]) # Q * M[0] * Q^t = anti-id
R = diagonalize_full_alternating_matrix(F,M[1]) # R * M[1] * R^t = anti-id
P = R.inverse()*Q
print("Check: P*M[0]*P^t == M[1] ", P*M[0]*P.transpose() == M[1])
print(P)

# Compute T (m*m) such that T*L_0 = L_1*P with L_i of size m *n and P size n*n
# Loop over all possible pairs (i, j) where i < j
for i in range(len(L_submatrix[0]) - 1):
    for j in range(i + 1, len(L_submatrix[0])):  # Ensure j > i
        try:
            # Get the i-th and j-th submatrices from L_submatrix[0]
            L0_i = L_submatrix[0][i]
            L0_j = L_submatrix[0][j]
            
            # Compute the pseudoinverse of L0_i
            L0_i_pseudo_inv = L0_i.pseudoinverse()

            # Compute T using L0_j and the pseudoinverse of L0_i
            T = (L0_j * P) * L0_i_pseudo_inv

            # Check if the transformation T satisfies the condition
            if T * L0_i != L0_j * P:
                # If the condition does not hold, continue to the next (i, j)
                continue
            
        except Exception as e:
            # If any error occurs (e.g., singular matrix), skip to the next pair (i, j)
            print(f"Error occurred at pair (i={i}, j={j}): {e}. Skipping to the next pair.")
            continue  # Go to the next (i, j) pair

        # If the condition holds, output the result
        print(f"\nCheck T*L[0][{i}] == L[0][{j}]*P: True, rank T is full? {T.rank() == m}")
        
        # Store the submatrices L_0_0 and L_0_1 as per the condition
        L_0_0 = L0_i
        L_0_1 = L0_j

        # Output the results and the matrices
        print("Condition holds! Storing L_0_0 and L_0_1:")
        print(f"L_0_0 (from index {i}):\n{L_0_0}")
        print(f"L_0_1 (from index {j}):\n{L_0_1}")

        # Break the inner loop (j loop) since we found a valid pair
        break
    else:
        # If no valid (i, j) pair found in this iteration, continue with the next i
        continue
    # Break the outer loop (i loop) as well since we have found a solution
    break

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
L1 = B * L_0_0

# Compute L using the second formula: L = A * L_1 * P
L2 = A * L_0_1

# Display results
print("Check if it yields the same L ",L1 == L2)

# Check if x, y in span(L_1 rows) then x*M[0]*y = 0
# Get the row space of L1
row_space_L1 = L1.row_space()

# Perform the check for all pairs of x, y in the row space of L1
def check_condition_2(M0,M1, row_space):
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
is_valid = check_condition_2(M[0], M[1],row_space_L1)

# Output the result
if is_valid:
    print("The condition x*M[0]*y = x*M[1]*y = 0 DOES hold for all x, y in the span of the rows of L1.")
else:
    print("The condition x*M[0]*y = x*M[1]*y = 0 DOES NOT hold for all x, y in the span of the rows of L1.")

print("L1 is of rank ", L1.rank())

# Check
# print(O_I.transpose())
# Get the row space of O_I.transpose()
row_space_O_I = O_I.transpose().row_space()
is_valid = check_condition_2(M[0], M[1],row_space_O_I)
# Output the result
if is_valid:
    print("The condition x*M[0]*y = x*M[1]*y = 0 DOES hold for all x, y in the span of the rows of O_I.")
else:
    print("The condition x*M[0]*y = x*M[1]*y = 0 DOES NOT hold for all x, y in the span of the rows of O_I.")

# Overal check if we found an oil space
# Perform the check for all pairs of x, y 
def check_condition_1(M0, row_space):
    basis = row_space.basis()
    for x in basis:
        for y in basis:
            # Compute x * M0 * y^T
            result_1 = x * M0 * y.column()
            if result_1 == 0:
                return True  # If any result is not zero, return False
    return False  # If all results are zero, return True

all_conditions_hold = True

# Iterate over all matrices M[i] and check the condition
for i in tqdm(range(m), ncols = 100, desc = "Checking oil space found or not ... "):
    M_i = M[i]  # Access the i-th matrix M[i]
    
    # Assume L1 is a matrix whose row space we are considering
    row_space_L1 = L1.row_space()

    # Check the condition for the current M_i
    is_valid = check_condition_1(M_i, row_space_L1)

    # Output the result for each M[i]
    if not is_valid:
        all_conditions_hold = False
    
# Final summary after checking all M[i]
Oil_subspace = row_space_L1
if all_conditions_hold:
    print("The condition x*M[i]*y = 0 holds for all M[i].")
    # Print the basis and dimension of the subspace
    # print("Oil subspace L found with basis:", Oil_subspace.basis())
    print("Dimension of the oil subspace L found:", print(Oil_subspace))
else:
    print("The condition x*M[i]*y = 0 does not hold for at least one M[i].")

# Check UOV vanish
# print(Oil_subspace)
print(check_uov_vanishing(F,Oil_subspace,[M[0],M[1]]))