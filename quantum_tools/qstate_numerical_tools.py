from sympy.physics.quantum import TensorProduct
from sympy import sqrt, Matrix
import numpy as np

def create_basis_ket_bra():
    # Basis state |0> and its conjugate transpose <0|
    ket0 = Matrix([[1], [0]])
    bra0 = ket0.T
    # Basis state |1> and its conjugate transpose <1|
    ket1 = Matrix([[0], [1]])
    bra1 = ket1.T
    return ket0, bra0, ket1, bra1

def create_ghz_state(n):
    if n < 1:
        raise ValueError("Number of qubits must be at least 1.")
    
    # GHZ state is a superposition of |00...0> and |11...1> states
    size = 2**n  # Total number of states in the n-qubit system
    
    # Create |00...0> state: ket0
    ket0 = np.zeros(size)
    ket0[0] = 1  # |00...0> corresponds to the first basis state
    
    # Create |11...1> state: ket1
    ket1 = np.zeros(size)
    ket1[-1] = 1  # |11...1> corresponds to the last basis state
    
    # Convert numpy arrays to sympy matrices
    ket0 = Matrix(ket0)
    ket1 = Matrix(ket1)
    
    # GHZ state: (|00...0> + |11...1>) / sqrt(2)
    ket_ghz = (ket0 + ket1) / sqrt(2)
    bra_ghz = ket_ghz.T  # Conjugate transpose of GHZ state
    return ket_ghz, bra_ghz

def tensor_n_times(state, n):
    """
    Returns the tensor product of the given state with itself n times.
    
    Parameters:
    state (Matrix): The initial quantum state.
    n (int): The number of times to apply the tensor product.
    
    Returns:
    Matrix: The resulting state after applying the tensor product n times.
    """
    # Start with the initial state
    result_state = state
    
    # Apply the tensor product n-1 times
    for _ in range(1, n):
        result_state = TensorProduct(result_state, state)
    
    return result_state

def create_initial_noisy_states(ket0, bra0, ket1, bra1, q, noise_type):
    """
    Generates initial quantum states under specified noise types.
    
    Parameters:
    ket0 (Matrix): The ket |0> state.
    bra0 (Matrix): The bra <0| state.
    ket1 (Matrix): The ket |1> state.
    bra1 (Matrix): The bra <1| state.
    q (Symbol): The quantum noise parameter.
    noise_type (str): The type of noise ('depolarizing', 'dephasing', 'bitflip').
    
    Returns:
    state1 (Matrix): First noisy state.
    state2 (Matrix): Second noisy state.
    state3 (Matrix): Third noisy state.
    state4 (Matrix): Fourth noisy state.
    
    Raises:
    ValueError: If an unsupported noise type is specified.
    """
    if noise_type == 'depolarizing':
        # Depolarizing noise affects both basis states and cross terms.
        state1 = (1 - q/2) * TensorProduct(ket0, bra0) + (q/2) * TensorProduct(ket1, bra1)
        state2 = (1 - q/2) * TensorProduct(ket1, bra1) + (q/2) * TensorProduct(ket0, bra0)
        state3 = (1 - q) * TensorProduct(ket0, bra1)
        state4 = (1 - q) * TensorProduct(ket1, bra0)
    elif noise_type == 'dephasing':
        # Dephasing noise primarily affects the off-diagonal terms.
        state1 = TensorProduct(ket0, bra0)
        state2 = TensorProduct(ket1, bra1)
        state3 = (1 - 2*q) * TensorProduct(ket0, bra1)
        state4 = (1 - 2*q) * TensorProduct(ket1, bra0)
    elif noise_type == 'bitflip':
        # Bitflip noise includes probabilities for flipping between basis states.
        state1 = (1 - q) * TensorProduct(ket0, bra0) + q * TensorProduct(ket1, bra1)
        state2 = (1 - q) * TensorProduct(ket1, bra1) + q * TensorProduct(ket0, bra0)
        state3 = (1 - q) * TensorProduct(ket0, bra1) + q * TensorProduct(ket1, bra0)
        state4 = (1 - q) * TensorProduct(ket1, bra0) + q * TensorProduct(ket0, bra1)
    else:
        raise ValueError("Unsupported noise type. Choose 'depolarizing', 'dephasing', or 'bitflip'.")
    
    return state1, state2, state3, state4

def compute_noisy_density_matrix(n, q, noise_type):
    """
    Computes the density matrix by taking the tensor product of the states for n qubits
    under the specified noise model and averaging them.
    
    Parameters:
    n (int): The number of qubits.
    q (Symbol): The quantum noise parameter.
    noise_type (str): The type of noise ('depolarizing', 'dephasing', 'bitflip').
    
    Returns:
    Matrix: The resulting averaged state matrix.
    """
    # Create basis states
    ket0, bra0, ket1, bra1 = create_basis_ket_bra()

    # Create initial noisy states
    state1, state2, state3, state4 = create_initial_noisy_states(ket0, bra0, ket1, bra1, q, noise_type)

    # Compute the tensor product of the state for n qubits and average them
    result_matrix = (tensor_n_times(state1, n) + tensor_n_times(state2, n) + 
                     tensor_n_times(state3, n) + tensor_n_times(state4, n)) / 2

    return result_matrix

def check_same_bits(number, n, n_anon):
    """
    Checks if the first n_anon bits of the binary representation of the number are the same. This function is to eliminate the state which collapsed when measurement in the generate anonymous entanglement protocol commence
    
    Parameters:
    number (int): The number to check.
    n (int): The total number of bits.
    n_anon (int): The number of anonymous bits to check.
    
    Returns:
    bool: True if all first n_anon bits are the same, False otherwise.
    """
    binary_string = f'{number:0{n}b}'
    
    # Check if the length of the binary string is less than n_anon
    if len(binary_string) < n_anon:
        return False
    
    # Slice the first n_anon bits
    first_n_bits = binary_string[:n_anon]
    
    # Check if all characters in the slice are the same as the first one
    return all(bit == first_n_bits[0] for bit in first_n_bits)

def fidelity_anonymous(result_matrix, n, n_anon):
    """
    Calculates the fidelity of the generated anonymous GHZ state.
    
    Parameters:
    result_matrix (Matrix): The result matrix from the noisy state computation.
    n (int): The number of qubits.
    n_anon (int): The number of anonymous bits.
    
    Returns:
    float: The calculated fidelity.
    """
    formula = 0
    for i in range(2**n):
        if check_same_bits(i, n, n_anon):
            formula += result_matrix[i, i]
    formula += result_matrix[0, 2**n-1] + result_matrix[2**n-1, 0]
    return formula / 2