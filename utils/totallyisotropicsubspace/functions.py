# Import all necessary Sage functionality
from sage.all import *

# Define a function anti identity matrix that has rank n = 2k
def create_anti_identity_matrix(Fq,k):
    # Define the 2x2 block
    A = matrix(Fq,[[0, 1], [1, 0]])

    # Create a list of n blocks, each of size 2x2
    blocks = [A for i in range(k)]

    anti_identity = block_diagonal_matrix(blocks)

    # Return the result
    return anti_identity


