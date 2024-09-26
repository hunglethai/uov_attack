# Import all necessary Sage functionality
from sage.all import *
from tqdm import tqdm
import itertools

# Define a function anti identity matrix that has rank n = 2k
def create_anti_identity_matrix(F: FiniteField,k: int) -> Matrix:
    """
    Generate a matrix like this
    
    [0 1 | 0 0]
    [1 0 | 0 0]
    -----------
    [0 0 | 0 1]
    [0 0 | 1 0]

    Args:
        F (FiniteField): Finite Field
        k (int): 1/2 dimension of the generated matrix 

    Returns:
        Matrix:     
                [0 1 | 0 0]
                [1 0 | 0 0]
                -----------
                [0 0 | 0 1]
                [0 0 | 1 0]
    """
    # Define the 2x2 block
    A = matrix(F,[[0, 1], [1, 0]])

    # Create a list of n blocks, each of size 2x2
    blocks = [A for i in range(k)]

    anti_identity = block_diagonal_matrix(blocks)

    # Return the result
    return anti_identity

# Function to check if a subspace is totally isotropic
def is_totally_isotropic(subspace_basis: list['Vector'], A: Matrix) -> bool:
    """
    Check if a basis of a vector subspace is totally isotropic

    Parameters:
    -----------
    subspace_basis: a list of vectors
    
    A: Matrix
    -----------

    Returns:
        bool: True/False
    """
    for v in subspace_basis:
        for w in subspace_basis:
            if (v * A * w) != 0:
                return False
    return True

# Brute-force search for totally isotropic subspace
def find_totally_isotropic_subspace(A: Matrix, F: FiniteField, dim_subspace: int) -> list['Vector']:
    """
    Finds a totally isotropic subspace of a given dimension for a bilinear form matrix A over a finite field F.
    
    Parameters:
    -----------
    A : sage.matrix.matrix_space.Matrix
        A square Sage matrix representing the bilinear or quadratic form. It should be symmetric and have dimensions n x n.
    
    F : sage.rings.finite_rings.finite_field.FiniteField
        The finite field over which the matrix A is defined.
    
    dim_subspace : int
        The dimension of the desired isotropic subspace. This should be less than or equal to the dimension of the matrix A (i.e., dim_subspace <= n).
    
    Returns:
    --------
    list of sage.modules.free_module_element.Vector:
        A list of vectors that form the basis for the totally isotropic subspace of the given dimension.
    
    Raises:
    -------
    ValueError:
        If the desired dimension is larger than the dimension of the matrix A or if other invalid conditions are met.
    """
    
    n = A.nrows()  # Dimension of the vector space
    
    # Generate all possible subspaces of given dimension
    vector_space = VectorSpace(F, n)
    all_combinations = list(vector_space.subspaces(dim_subspace))
    # print(all_combinations)
    
    for subspace in tqdm(all_combinations, desc="Searching for totally isotropic subspace"):
        # Get a basis for the subspace
        subspace_basis = subspace.basis()
        
        # Check if this subspace is totally isotropic
        if is_totally_isotropic(subspace_basis, A):
            return subspace_basis  # Return the basis of the isotropic subspace
    
    return None  # No totally isotropic subspace found

# Function to generate matrix from basis vectors using a loop
def generate_matrix_from_basis(F: FiniteField,basis: list['Vector']) -> Matrix:
    """
    Generate a matrix from a list of vectors in a basis

    Args:
        F (_type_): Finite field
        basis (_type_): list of vectors

    Returns:
        representation_matrix: A matrix that represents the basis
    """
    # Create an empty matrix with the same number of rows as the vectors and as many columns as there are basis vectors
    rows = len(basis[0])  # The length of each vector gives the number of rows
    cols = len(basis)     # The number of basis vectors gives the number of columns
    representation_matrix = matrix(F, rows, cols)  # Initialize the matrix over F
    
    # Populate the matrix by adding each vector as a column
    for j, vec in enumerate(basis):
        for i in range(rows):
            representation_matrix[i, j] = vec[i]  # Place each element in the appropriate position
    
    return representation_matrix

# Function to check if isotropic_matrix.transpose() * L.transpose() * A * L * isotropic_matrix == 0
def check_condition(isotropic_matrix: Matrix, A: Matrix, L: Matrix,F: FiniteField):
    """
    Check if isotropic_matrix.transpose() * L.transpose() * A * L * isotropic_matrix == 0

    Args:
        isotropic_matrix (Matrix): Target matrix to check
        A (Matrix): bilinear form or quadratic form representation matrix
        L (Matrix): A change of basis matrix if any
        F (FiniteField): Finite field

    Returns:
        True/False: True if == 0
    """
    # Ensure L is not a zero matrix
    if L.is_zero():
        return False  # L is zero, reject it
    # Proceed with the original condition check
    n, k = isotropic_matrix.nrows(), isotropic_matrix.ncols()
    result = isotropic_matrix.transpose() * L.transpose() * A * L * isotropic_matrix
    return (result.is_zero())

# Brute-force search for matrix L such that isotropic_matrix.transpose() * L.transpose() * A * L * isotropic_matrix == 0
def find_L_for_condition(isotropic_matrix: Matrix, A: Matrix, F: FiniteField) -> list[Matrix]:
    """
    Brute-force search for matrix L such that isotropic_matrix.transpose() * L.transpose() * A * L * isotropic_matrix == 0

    Args:
        isotropic_matrix (Matrix): A matrix
        A (Matrix): bilinear form or quadratic form representation matrix
        F (FiniteField): Finite field

    Returns:
        list[Matrix]: List of L
    """
    n = A.nrows()  # Dimension n
    q = F.order()

    # Generate all possible matrices of size n * n over F
    # This will be slow for large n and q, as we're brute-forcing it
    num_matrices = q ** (n * n)
    print(f"Brute-forcing {num_matrices} matrices of size {n}x{n} over GF({q})...")
    
    valid_L_matrices = []  # List to store all valid L matrices found

    # Progress bar using tqdm
    for entries in tqdm(cartesian_product_iterator([F]*n*n), total=int(num_matrices), desc="Searching for L"):
        # Reshape the list of entries into a matrix L of size n * n
        L = matrix(F, n, n, entries)
        
        # Check if L satisfies the condition
        if check_condition(isotropic_matrix, A, L, F):
            valid_L_matrices.append(L)  # Add valid L to the list

    # Return the list of valid L matrices
    return valid_L_matrices

# Brute-force search for matrices L such that L.transpose() * A * L == 0
def brute_force_search_isotropic_matrices(A: Matrix, F: FiniteField, m: int) -> Matrix:
    """
    Brute-force search for a matrix L of size n x m such that L^T * A * L == 0.

    Parameters:
    -----------
    A : Matrix
        A square matrix of size n x n representing the bilinear form or quadratic form.
    F : FiniteField
        The finite field over which the entries of matrix A, L are defined.
    m : int
        The number of columns of matrix L.

    Returns:
    --------
    Matrix or None:
        Returns a matrix L of size n x m that satisfies L^T * A * L == 0,
        or None if no such matrix is found.
    """
    n = A.nrows()
    valid_matrices = []
    
    # Total combinations
    total_combinations = len(F) ** (n * m)

    # Iterate over all possible values for matrix L in the finite field
    for values in tqdm(itertools.product(F, repeat=n * m), total=total_combinations, desc="Searching for valid matrices"):
        # Create a matrix L from the current combination of values
        L = matrix(F, n, m, values)
        
        # Check if L^T * A * L == 0
        if (L.transpose() * A * L).is_zero():
            valid_matrices.append(L)  # Add the valid matrix to the list

    return valid_matrices  # Return the list of valid matrices

# Check if a matrix A is an alternating matrix
def check_alternating_matrix(A: Matrix, F: FiniteField) -> bool:
    """
    Check if a matrix is an alternating matrix

    Args:
        A (Matrix): Any matrix
        F (FiniteField): Finite field

    Returns:
        bool: True/False
    """
    n_rows = A.nrows()
    n_cols = A.ncols()
    
    # The matrix must be square
    if n_rows != n_cols:
        return False
    
    # Check if all elements in the diagonal is zero or not:
    for i in range(n_rows):
        if A[i,i] != 0:
            return False
        
    # Check if the matrix is skew-symmetric, i.e. A = - A.transpose()
    for i in range(n_rows):
        for j in range(i+1,n_cols):
            if A[i,j] != A[j,i]:
                return False
            
    return True

# Diagonalize top-left 2x2 block
def diagonalize_2x2_alternating_matrix(A: Matrix, F: FiniteField) -> Matrix:
    """
    Diagonalize the top-left 2x2 block

    Args:
        A (Matrix): Matrix
        F (FiniteFiel): Finite field

    Returns:
        Matrix: Diagonalized-top-left-block matrix
    """
    # Check if the input matrix is a full-rank alternating matrix
    if check_alternating_matrix(A,F) == False:
        print("A is not a full-rank alternating matrix")
        return None    
    
    if A.is_zero():
        return matrix.identity(F,A.nrows())
    
    # Number of rows and columns
    n_rows = A.nrows()
    n_cols = A.ncols()
    
    # Check of A[0,1] != 0
    if A[0,1] != 0:   
        # Create an identity matrix of size (n_rows - 2) x (n_cols - 2)
        I = matrix.identity(F,n_rows -2)
        
        # Create top-left matrix
        L_1 = matrix(F,[[1,0],[0,A[0,1].inverse()]])
        
        # Create top-right matrix
        L_2 = zero_matrix(F,2,n_cols-2)
        
        # Create bottom-left matrix
        list_values = []
        for i in range(2,n_rows):
            list_values.append(A[0,i]*A[0,1].inverse())
            
        L_3 = matrix(F,[[0]*(n_rows - 2),list_values])
        L_3 = L_3.transpose()
        
        # Construct L
        L = block_matrix(F,[[L_1,L_2],[L_3,I]],subdivise = False)
        
        return L
    
    if A[0,1] == 0:
        print("\n Switching ... \n")
        # Search for any A[0,i] != 0, it must exist at least one since A is a full-rank matrix
        for i in range(0,n_cols):
            if A[0,i] != 0:
                switch_index = i # This is the switching index
                break

        # Create a column switching matrix of size similar to A
        I = matrix.identity(F,n_rows)
        
        # Swap columns i and j
        I.swap_columns(i, 1)
        
        # Return the matrix A_new that has column 2 and i+1 switched
        A_switched = I*A*I.transpose()
        
        # Create bottom-righ matrix i.e. an identity matrix of size (n_rows - 2) x (n_cols - 2)
        K = matrix.identity(F,n_rows -2)
        
        # Create top-left matrix
        L_1 = matrix(F,[[1,0],[0,A_switched[0,1].inverse()]])
        
        # Create top-right matrix
        L_2 = zero_matrix(F,2,n_cols-2)
        
        # Create bottom-left matrix
        list_values = []
        for i in range(2,n_rows):
            list_values.append(A_switched[0,i]*A_switched[0,1].inverse())
            
        L_3 = matrix(F,[[0]*(n_rows - 2),list_values])
        L_3 = L_3.transpose()
        
        # Construct L
        L = block_matrix(F,[[L_1,L_2],[L_3,K]],subdivise = False)
        
        return L*I # Return switched L
    
# Diagonalize an alternating matrix
def diagonalize_full_alternating_matrix(A: Matrix, F: FiniteField) -> Matrix:
    """
    Return a matrix L such that L is invertible and L.transpose()*A*L = 
    
    [0 1 0 0]
    [1 0 0 0]
    [0 0 0 1]
    [0 0 1 0]

    Args:
        A (Matrix): an alternating matrix i.e. A is skew-symmetric and has zeros in its diagonal
        F (FiniteField): Finite Field

    Returns:
        Matrix: Matrix L
    """
    
    # Rows and Columns
    n_rows = A.nrows()
    n_cols = A.ncols()
    # A_original = A.copy()
    
    # Diagonalizing block by block
    list_i_A_L = []
    
    for i in tqdm(range(n_rows),ncols = 100,desc="Diagonalizing block by block ..."):
        if i == 0:
            # Get the diagonalizing matrix of the first block
            L = diagonalize_2x2_alternating_matrix(A,F)
            
            # Get new A
            A = L*A*L.transpose()
            
            # Extract a submatrix of size (n_rows - 1) x (n_cols - 1)
            A = A.submatrix(1, 1, A.nrows()-1, A.ncols()-1)
            
            # Append to the list
            list_i_A_L.append((i,A,L))
            
        if i != 0:
            A = list_i_A_L[i-1][1]

            # Get the diagonalizing matrix of the first block as bottomright corner
            L = diagonalize_2x2_alternating_matrix(A,F)
            
            # Get new A
            A = L*A*L.transpose()
            
            # Get the identity matrix of size i * i as topleft corner
            I = matrix.identity(F,i)

            # Construct the block matrix
            L = block_diagonal_matrix(I,L)

            # Extract a submatrix of size (n_rows - 1) x (n_cols - 1)
            A = A.submatrix(1, 1, A.nrows()-1, A.ncols()-1)
                        
            # Append to the list
            list_i_A_L.append((i,A,L))
                        
    # Initialize the result_L as an identity matrix of the appropriate size
    result_L = matrix.identity(F, n_rows)   
    
    # Iterate through list_i_A_L and multiply all matrices L
    for _, _, L in list_i_A_L[::-1]:
        result_L *= L  # Multiply result by each matrix L

    # The result_L now contains the product of all matrices in list_L such that result_L*A*result_L.transpose() is almost anti identity
    
    # Final permutation matrix to get anti identity matrix
    list_P = []
    
    for i in range(2, n_rows, 2):
        P = identity_matrix(n_rows)
        P[i, i-2] = 1   # Add a 1 in the row two steps above
        list_P.append((i,P))
    
    result_P = matrix.identity(F, n_rows) 
    for _, P in list_P[::-1]:
        result_P *= P 

    return result_P*result_L
