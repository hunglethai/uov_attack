# Import all necessary Sage functionality
from sage.all import *
from tqdm import tqdm

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
    A = matrix(Fq,[[0, 1], [1, 0]])

    # Create a list of n blocks, each of size 2x2
    blocks = [A for i in range(k)]

    anti_identity = block_diagonal_matrix(blocks)

    # Return the result
    return anti_identity

# Function to check if a subspace is totally isotropic
def is_totally_isotropic(subspace_basis: list[Vector], A: Matrix) -> bool:
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
def find_totally_isotropic_subspace(A: Matrix, F: FiniteField, dim_subspace: int) -> list[Vector]:
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
def generate_matrix_from_basis(F: FiniteField,basis: list[Vector]) -> Matrix:
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


