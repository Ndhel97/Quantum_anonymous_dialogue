import numpy as np
import cmath
import random
import math

def generate_parties(n):
    """
    Generate a dictionary of party identifiers for qubits assignment.

    Parameters:
    n (int): A positive integer indicating number of participants in the network.
             Must be at least 2.

    Returns:
    A tuple containing:
        - dict: A dictionary mapping party identifiers to their corresponding integer values.
        - int: The total number of qubits.

    Raises:
    ValueError: If n is less than 2.
    """
    if n < 2:
        raise ValueError("n must be at least 2.")
    
    p = {
        'a0': 0,        'a1': 1, 
        'b0': 2,        'b1': 3,
        'a2': 2 * n,    'a3': 2 * n + 1,
        'b2': 2 * n + 2,'b3': 2 * n + 3,
        'A1': 4 * n,    'A2': 4 * n + 1,
        'B1': 4 * n + 2,'B2': 4 * n + 3
    }

    number_of_qubits = 4 * n + 4
    
    return p, number_of_qubits


def generate_probabilities(sample_probs):
    """
    Generate an array of probabilities (noise parameter) ranging from 0 to 1, inclusive.

    Parameters:
    sample_probs (int): The number of sample probabilities to generate.

    Returns:
    np.ndarray: An array of probabilities ranging from 0 to 1, inclusive.
    """
    if sample_probs <= 0:
        raise ValueError("sample_probs must be a positive integer.")
    
    return np.linspace(0, 1, num=sample_probs)

def generate_uni_sample(n):
    """
    Generate a sample of n pure states with normalized probability amplitudes.
    The coefficients alpha_0 and alpha_1 of the quantum states are generated uniformly on the unit sphere.

    Parameters:
    n (int): The number of samples to generate.

    Returns:
    np.ndarray: A numpy array of shape (n, 2) containing the generated coefficients.
    """
    uni_sample = np.ndarray(shape=(n, 2), dtype=complex)
    for i in range(n):
        gamma_1 = random.uniform(0.0, 1.0)
        gamma_2 = random.uniform(0.0, 1.0)
        theta = math.acos(1 - 2 * gamma_1)
        phi = 2 * math.pi * gamma_2
        
        alpha_0 = math.cos(theta / 2)
        alpha_1 = cmath.exp(1j * phi) * math.sin(theta / 2)
        
        vector = np.array([alpha_0, 0, 0, alpha_1])
        magnitude = np.linalg.norm(vector)
        
        uni_sample[i][0] = alpha_0 / magnitude
        uni_sample[i][1] = alpha_1 / magnitude
           
    return uni_sample

def clean_complex_array(arr, threshold=1e-10):
    """
    Clean an array of complex numbers: round near-zero reals and imaginaries to zero,
    and format numbers to avoid unnecessary scientific notation.

    Parameters:
    arr (np.ndarray): Numpy array of complex numbers.
    threshold (float): Magnitude below which numbers are rounded to zero. Default is 1e-10.

    Returns:
    np.ndarray: A cleaned numpy array with the same shape as `arr`.
    """
    def clean_number(x):
        real = 0.0 if abs(x.real) < threshold else np.round(x.real, 10)
        imag = 0.0 if abs(x.imag) < threshold else np.round(x.imag, 10)
        return complex(real, imag)
    
    # Vectorize the clean_number function
    clean_func = np.vectorize(clean_number)
    
    # Apply the cleaning function to each element in the array
    clean_arr = clean_func(arr)
    
    return clean_arr

def form_dm(alpha_0, alpha_1):
    """
    Form the density matrix based on the coefficients of the quantum state.

    Parameters:
    alpha_0 (complex): Coefficient for the |0⟩ state.
    alpha_1 (complex): Coefficient for the |1⟩ state.

    Returns:
    np.ndarray: The cleaned density matrix.
    """
    vector = np.array([alpha_0, 0, 0, alpha_1])
    rho = np.outer(vector, np.conj(vector.T))
    dm = clean_complex_array(rho)
    return dm

def generate_ghz_dm(n):
    """
    Generate the density matrix for an n-qubit GHZ state.

    Parameters:
    n (int): Number of qubits.

    Returns:
    np.ndarray: The density matrix of the GHZ state.
    """
    # Calculate the dimension of the matrix, which is 2^n
    dim = 2**n
    
    # Create an empty matrix of complex type
    dm_ghz = np.zeros((dim, dim), dtype=complex)
    
    # Set the amplitude for the GHZ state components
    amplitude = 1 / np.sqrt(2)
    
    # Set the first and last elements of the matrix
    dm_ghz[0, 0] = amplitude * np.conj(amplitude)
    dm_ghz[0, -1] = amplitude * np.conj(amplitude)
    dm_ghz[-1, 0] = amplitude * np.conj(amplitude)
    dm_ghz[-1, -1] = amplitude * np.conj(amplitude)
    
    return dm_ghz

def update_party_indices(p):
    """
    Update the dictionary of party identifiers to new indices. Add qubits for the messages.

    Parameters:
    p (dict): The dictionary of party identifiers to be updated.

    Updates the provided dictionary with new mappings for the party identifiers.
    """
    p.update({
        'a1': 0, 'b0': 1, 'b1': 2,
        'a2': 3, 'a3': 4, 'b2': 5,
        'A1': 6, 'A2': 7, 'B1': 8, 'B2': 9
    })

