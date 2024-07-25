#TODO:  Add validation functions for matrix layers / .X / .raw 
import numpy as np
from scipy.sparse import issparse


def is_integers(array_like):
    """
    Purpose: 
    Checks if a given array contains only integer values.
    
    Parameters:
        array_like: The array to be checked.
    
    Returns: 
        True if all elements are integers, otherwise False.
    """

    elements_to_check = 100000
    if issparse(array_like):
        data = array_like.data[:elements_to_check]
    else:
        rows_to_check = np.ceil(elements_to_check/array_like.shape[1]).astype('int')
        data = array_like[:rows_to_check,:]
    result = np.all(np.mod(data, 1) == 0)
        
    return result
