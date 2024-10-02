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

def is_skew_symmetric(F,M: Matrix) -> bool:
    # Get the number of rows (and columns, assuming it's square)
    n = M.nrows()
    
    # Iterate through the matrix and check the skew-symmetric condition
    for i in range(n):
        for j in range(i, n):
            if M[i, j] != -M[j, i]:
                return False
    return True

# Generate a matrix that vanishes on the oil space
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

# Generate a list of matrices where each matrix is the sum of the corresponding matrix in the input list and its transpose
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

# Create a function that checks if a subspace L vanishes in a list of matrices M[i].
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

        
    
    