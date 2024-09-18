# Import all necessary Sage functionality
from sage.all import *
from tqdm import tqdm

# Define a function anti identity matrix that has rank n = 2k
def create_anti_identity_matrix(Fq,k):
    # Define the 2x2 block
    A = matrix(Fq,[[0, 1], [1, 0]])

    # Create a list of n blocks, each of size 2x2
    blocks = [A for i in range(k)]

    anti_identity = block_diagonal_matrix(blocks)

    # Return the result
    return anti_identity

# Function to check if a subspace is totally isotropic
def is_totally_isotropic(subspace_basis, A):
    for v in subspace_basis:
        for w in subspace_basis:
            if (v * A * w) != 0:
                return False
    return True

# Brute-force search for totally isotropic subspace
def find_totally_isotropic_subspace(A, F, dim_subspace):
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
def generate_matrix_from_basis(F,basis):
    # Create an empty matrix with the same number of rows as the vectors and as many columns as there are basis vectors
    rows = len(basis[0])  # The length of each vector gives the number of rows
    cols = len(basis)     # The number of basis vectors gives the number of columns
    isotropic_matrix = matrix(F, rows, cols)  # Initialize the matrix over F
    
    # Populate the matrix by adding each vector as a column
    for j, vec in enumerate(basis):
        for i in range(rows):
            isotropic_matrix[i, j] = vec[i]  # Place each element in the appropriate position
    
    return isotropic_matrix

# Function to check if isotropic_matrix.transpose() * L.transpose() * A * L * isotropic_matrix == 0
def check_condition(isotropic_matrix, A, L,F):
    # Ensure L is not a zero matrix
    if L.is_zero():
        return False  # L is zero, reject it
    # Proceed with the original condition check
    n, k = isotropic_matrix.nrows(), isotropic_matrix.ncols()
    result = isotropic_matrix.transpose() * L.transpose() * A * L * isotropic_matrix
    return (result.is_zero())

# Brute-force search for matrix L
def find_L_for_condition(isotropic_matrix, A, F):
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


