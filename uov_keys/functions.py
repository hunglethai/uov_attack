from sage.all import *
from tqdm import tqdm

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
    
    for i in tqdm(range(m), ncols =  100, desc = "Generating public keys ... "):
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
        
        
    
    