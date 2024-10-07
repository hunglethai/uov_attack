from sage.all import *
from tqdm import tqdm
from itertools import combinations

# Generate a random full rank matrix
def random_full_rank_matrix(F: FiniteField, m: int, n: int):
    """
    Generate a random matrix of size m x n which is full rank
    """
    while True:
        # Generate a random m x n matrix over the specified field (default is rational field QQ)
        A = random_matrix(F, m, n)
        
        # Check if the rank is full (i.e., min(m, n))
        if A.rank() == min(m, n):
            return A

def random_upper_triangular_matrix(F: FiniteField, m: int, n: int) -> Matrix:
    """
    Generates a random upper-triangular matrix of size m x n over the specified field.
    
    m: number of rows
    n: number of columns
    F: Finite field
    
    Returns:
        A random upper-triangular matrix of size m x n.
    """
    # Initialize an empty matrix of size m x n
    M = Matrix(F, m, n)
    
    # Fill in the upper triangular part with random elements
    for i in range(m):
        for j in range(i, n):
            M[i, j] = F.random_element()  # Fill the upper triangular part with random values
    
    # Elements below the diagonal are automatically zero, as we do not modify them
    
    return M

# Generate oil subspace matrix i.e. the secret key
def generate_oil_subspace_matrix(F: FiniteField, m: int, n: int) -> Matrix:
    """
    F: Finite Field
    m: oil variables i.e o
    n: oil and vinegar variables i.e. o + v
    
    Return 
    A matrix which is the basis of the oil subspace
    """
    A = random_full_rank_matrix(F,m,n-m)
    I = matrix.identity(F,m)
    O = block_matrix(1,2,[A,I], subdivise = False)
    
    return A.transpose(), O.transpose()

# "Upper()" function as in UOV specs
def upper(F: FiniteField, M: Matrix) -> Matrix:
    """
    Returns the unique upper triangular matrix M' such that M' - M is skew-symmetric.
    
    M: A square matrix
    """
    # Step 1: Calculate M + M^T (transpose of M)
    M_plus_MT = M + M.transpose()
    
    # Step 2: Create a matrix X with the lower triangular part of M + M^T
    n = M.nrows()
    X = Matrix(F, n, n)
    
    for i in range(n):
        for j in range(i+1):
            # Copy lower triangular part (including diagonal)
            X[i, j] = M_plus_MT[i, j]
    
    # Return the matrix X
    return X

# Check if a matrix is indeed skew-symmetric
def is_skew_symmetric(F,M: Matrix) -> bool:
    # Get the number of rows (and columns, assuming it's square)
    n = M.nrows()
    
    # Iterate through the matrix and check the skew-symmetric condition
    for i in range(n):
        for j in range(i, n):
            if M[i, j] != -M[j, i]:
                return False
    return True

# Generate a matrix that vanishes on the oil space i.e. public key matrices
def generate_public_matrices(F: FiniteField, m: int, n: int, O: Matrix) -> list:
    """
    F: Finite Field
    m: oil variables i.e o
    n: oil and vinegar variables i.e. o + v
    O: oil subspace
    
    Return 
    A list of matrices which vanish on O
    """
    list_P = []
    
    for i in range(m):
        P_1 = random_upper_triangular_matrix(F,n-m,n-m)
        P_2 = random_matrix(F,n-m,m)
        P_3 = upper(F,-O.transpose()*P_1*O - O.transpose()*P_2)
        P_4 = zero_matrix(F,m,n-m)
        P = block_matrix(2,2,[[P_1,P_2],[P_4,P_3]])
        list_P.append(P)
    
    return list_P

# Generate a list of matrices where each matrix is the sum of the corresponding matrix in the input list and its transpose i.e M = P + P^t i.e symplectic forms UOV
def generate_list_M(list_P):
    """
    Generate a list of matrices where each matrix is the sum of the corresponding 
    matrix in the input list and its transpose.
    
    Parameters:
    ----------
    list_P : list of matrices

    Returns:
    -------
    list_M : list of matrices
        A list of matrices where each matrix M is the sum of the corresponding 
        matrix P in `list_P` and its transpose P^T.
    """
    
    # Create an empty list to store the matrices in list_M
    list_M = []
    
    # Iterate over each matrix in list_P
    for P in list_P:
        # Compute P + P^T (transpose of P)
        M = P + P.transpose()
        # Append the resulting matrix to list_M
        list_M.append(M)
    
    return list_M

# Check if a subspace A (having basis X) is invariant under the linear transformation T
def is_invariant_subspace(F, X, T):
    """
    Check if the subspace A spanned by X is invariant under the linear transformation T.
    
    Args:
        F: The finite field over which the matrices are defined.
        X: A matrix of size m x n representing the basis of the subspace A (with m basis vectors).
        T: A matrix of size n x n representing the linear transformation.

    Returns:
        True if the subspace A is invariant under T, False otherwise.
    """
    # Create the subspace A spanned by the rows of X
    A = span(X.rows())

    # Iterate over each basis vector in X
    for v in X.rows():
        # Apply the transformation T to the vector v
        T_v = v*T

        # Check if the transformed vector T_v lies in the span of the basis of A
        if not T_v in A:
            return False

    # If all transformed vectors lie in A, the subspace is invariant under T
    return True

# Create a function that checks if a subspace L vanishes in a list of matrices M[i] i.e. UOV public keys
def check_uov_vanishing(F: FiniteField, L , M_list: list) -> bool:
    """
    Check if for all row vectors x, y in subspace L and for all matrices M in M_list,
    we have x * M * y^T = 0.
    
    Args:
        F: A finite field.
        L: A subspace of vectors over F (as a MatrixSpace or a list of row vectors).
        M_list: A list of matrices over F.

    Returns:
        True if for all x, y in L, x * M * y^T = 0 for each M in M_list, False otherwise.
    """
    # Get the basis matrix of L 
    basis_matrix = L.basis_matrix()

    # Iterate over all pairs of row vectors (x, y) from the basis matrix
    for M in tqdm(M_list, ncols = 100, desc = "Check if UOV public keys vanish on Oil subpace ... "):
        # Check that for all row vectors x, y in the row space, x * M * y^T = 0
        for i in range(basis_matrix.nrows()):
            for j in range(basis_matrix.nrows()):
                x = basis_matrix.row(i)  # row vector x
                y = basis_matrix.row(j)  # row vector y
                
                # Calculate the result of x * M * y^T
                result = x * M * y.column()  # y.column() converts y to a column vector
                
                # Ensure the result is treated as a scalar and compare with 0
                if result != 0:
                    return False
    return True

# Visualize matrices where non-zero elements are represented by black dots and zero elements by white dots
def color_matrix(matrix: Matrix):
    """
    Function to print a large matrix in a compact way, using:
    - Filled squares (■) for non-zero elements.
    - Empty squares (□) for zero elements.
    
    Args:
    matrix (Matrix): A matrix over a finite field.

    Prints:
    A compact grid-based representation of the matrix.
    """
    for i in range(matrix.nrows()):
        row_str = ""
        for j in range(matrix.ncols()):
            if matrix[i, j] == 0:
                row_str += "□ "  # Empty square for zero
            else:
                row_str += "■ "  # Filled square for non-zero
    print(row_str.strip())  # Print each row of squares

# Construct the isotropic subspace basis from a symplectic basis noted that if m = n//2 it is a Lagrangian        
def compute_isotropic_subspace_basis(F: FiniteField, L: list, M: list, m: int, n: int, lag_or_not: bool, o: int = None):
    """
    This function computes o-dimensional isotropic subspace basis form a symplectic basis.

    Args:
        L (list): A list of matrices which are the symplectic basis, from which submatrices will be extracted.
        M (list): A list of matrices that are symplectic forms.
        n (int): Dimension of the Finite Field i.e. n = 2o
        m (int): Number of rows in the isotropic subspace basis i.e. the isotropic subspace is m-dimensional. If m = n//2 we are dealing with lagrangians.
        F: A finite field or identity matrix to check against.
        lag_or_not: We are computing o-dimensional isotropic subspace basis or not ? i.e we are computing Lagrangians or not ?
            If True, the function will automatically set o = n // 2.
            If False, it will ask the user to input the value for o.
        o (int, optional): The number of rows to select. Will be automatically set if lag_or_not is True.

    Returns:
        matrices_list: A list containing all o-dimensional isotropic subspace basis that are full-rank of m symplectic from of dimension n.
        all_check_pass (bool): Boolean flag indicating if all submatrices pass the Lagrangian check.
    """
    # Set o based on lag_or_not flag
    if lag_or_not:
        o = n // 2
        print(f"We are computing Lagrangian basis, setting o to {o} (n // 2)")
    elif o is None:  # If lag_or_not is False and no value for o is provided
        o = int(input("We are computing isotropic basis of dimension o, please enter the value of o: "))

    # Initialize a list to store all submatrices of size o x n that meet the criteria
    L_submatrix = []
    
    # Iterate over all possible matrices in L
    for l in tqdm(range(m), ncols=100, desc=f"Computing all {o}-dimensional isotropic subspaces of each of {m} symplectic basis ... "):
        L_list = []
        
        # Loop over possible k values
        for k in range(o + 1):
            first_half_rows = range(n//2)
            second_half_rows = range(n//2, n)
            
            # Choose k rows from the first half
            for first_set in combinations(first_half_rows, k):
                # Exclude the corresponding second half rows
                forbidden_rows = [u + n//2 for u in first_set]
                valid_second_half_rows = [row for row in second_half_rows if row not in forbidden_rows]
                
                # Choose o - k rows from the second half
                for second_set in combinations(valid_second_half_rows, o - k):
                    # Combine selected rows
                    selected_rows = list(first_set) + list(second_set)
                    
                    # Extract submatrix
                    submatrix = L[l][selected_rows, :]
                    
                    # Check for full rank using Sage's rank function
                    if submatrix.rank() == o:
                        L_list.append(submatrix)
        
        L_submatrix.append(L_list)
    
    # Check if the submatrices are isotropic and full rank
    all_check_pass = True
    for i in range(m):
        for j in range(len(L_submatrix[i])):
            if (L_submatrix[i][j] * M[i] * L_submatrix[i][j].transpose()).is_zero() == False or (L_submatrix[i][j].rank() != L_submatrix[i][j].nrows()):
                all_check_pass = False
                break  # Exit loop if any check fails
    
    # Only return and print results if all checks pass
    if all_check_pass:
        print(f"All submatrices are {o}-dimensional isotropic subspace basis and full rank.")
        return L_submatrix
    else:
        print("Some submatrices failed the check. No output returned.")
        return None

# Computes the matrix T such that T * L_0 = L_1 * P for two submatrices L_0 and L_1 are in matrices_list
def compute_transformation_T(matrices_list: list, P: Matrix):
    """
    Computes the matrix T such that T * L_0 = L_1 * P for two submatrices L_0 and L_1 of size m * n, 
    where P is of size n * n. The function will find the pair (L_0, L_1) that satisfies the condition
    and return T and the submatrices if a solution is found.

    Args:
        matrices_list: A list of submatrices from which L_0 and L_1 are extracted.
        P: The matrix P of size n * n used in the transformation.

    Returns:
        T (matrix or None): The transformation matrix T if a valid one is found, otherwise None.
        L_0_0 (matrix or None): The submatrix L_0_0 if a solution is found, otherwise None.
        L_0_1 (matrix or None): The submatrix L_0_1 if a solution is found, otherwise None.
        found_solution (bool): Flag indicating if a solution was found.
    """
    
    # Initialize a flag to track if a solution is found
    found_solution = False
    T = None
    L_0_0 = None
    L_0_1 = None
    m = matrices_list[0].nrows()
    n = matrices_list[0].ncols()

    # Loop over all possible pairs (i, j)
    for i in tqdm(range(len(matrices_list)), ncols=100, desc="Computing T such that T*L_0 = L_1*P ..."):
        for j in range(len(matrices_list)):
            try:
                # Get the i-th and j-th submatrices from matrices_list
                L0_i = matrices_list[i]
                L0_j = matrices_list[j]
                L0_j_P = L0_j*P

                # Compute the last m columns of L0_i and L0_j_P
                L0_i_tail = L0_i[:, -m:]
                L0_j_P_tail = L0_j_P[:, -m:]

                # Compute T = L0_j_P_tail * L0_i_tail.inverse() 
                T = L0_j_P_tail * L0_i_tail.inverse() 

                # Check if the transformation T satisfies the condition
                if T * L0_i != L0_j * P:
                    # If the condition does not hold, continue to the next (i, j)
                    continue
                
            except Exception as e:
                # If any error occurs (e.g., singular matrix), skip to the next pair (i, j)
                # print(e)
                continue  # Go to the next (i, j) pair

            # If the condition holds, output the result
            print(f"\nCheck T*L[0][{i}] == L[0][{j}]*P: {T * L0_i == L0_j * P}, rank T is full? {T.rank() == T.nrows()}")
            
            # Store the submatrices L_0_0 and L_0_1 as per the condition
            L_0_0 = L0_i
            L_0_1 = L0_j

            # Output the results and the matrices
            print("Condition holds!")
            print(f"L_0_0 (from index {i})", L_0_0, "\n")
            print(f"L_0_1 (from index {j})", L_0_1, "\n")

            # Set the flag to True since a solution was found
            found_solution = True

            # Break the inner loop (j loop) since we found a valid pair
            break
        if found_solution:
            # Break the outer loop (i loop) as well since we have found a solution
            break

    # If no solution was found after all iterations, print "No solution"
    if not found_solution:
        print("No solution found !")
    else:
        print(f"Check T*L_0_0 == L_0_1*P ? {T * L_0_0 == L_0_1 * P}")

    # Return the results
    return T, L_0_0, L_0_1
