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

