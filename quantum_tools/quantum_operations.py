import netsquid as ns
from netsquid.qubits.qubitapi import *
from netsquid.qubits.operators import *
from netsquid.qubits.qformalism import QFormalism
import random


def apply_noise(qubit, noise, prob):
    """
    Apply noise to a given qubit based on the specified noise type and probability.

    Parameters:
    qubit (Qubit): The qubit to which noise will be applied.
    noise (str): The type of noise to apply. Supported types are "depolarize", "dephase", and "bitflip".
    prob (float): The probability of the noise affecting the qubit. Must be in the range [0, 1].

    Raises:
    ValueError: If an unsupported noise type is provided.

    Noise Types:
    - "depolarize": Applies depolarizing noise, which replaces the qubit's state with a completely mixed state
      with the given probability.
    - "dephase": Applies dephasing noise, which introduces phase damping (randomly flipping the phase) with
      the given probability.
    - "bitflip": Applies bit-flip noise, which flips the qubit's state (|0⟩ to |1⟩ and vice versa) with the given
      probability.
    """
    if noise == "depolarize":
        ns.qubits.depolarize(qubit, prob=prob)
    elif noise == "dephase":
        ns.qubits.dephase(qubit, prob=prob)
    elif noise == "bitflip":
        ns.qubits.apply_pauli_noise(qubit, (1-prob, prob, 0, 0))
    else:
        raise ValueError("Unsupported noise type")

def generate_anonymous(q, p, n, noise, GHZ_party, GHZ_3):
    """
    Generate the anonymous state using a specified noise model and calculate fidelities.

    Parameters:
    q (float): The noise probability parameter.
    p (dict): A dictionary mapping party identifiers to their corresponding qubit indices.
    n (int): Number of participants in the network.
    noise (str): The type of noise to apply. Supported types are "depolarize", "dephase", and "bitflip".
    GHZ_party (np.ndarray): The GHZ state for all party qubits.
    GHZ_3 (np.ndarray): The 3-qubit GHZ state.

    Returns:
    tuple: A tuple containing:
        - fidelity_anonymous (float): The fidelity of the anonymous state with the 3-qubit GHZ state.
        - fidelity_ghz (float): The fidelity of the generated GHZ state with the target GHZ state.
        - dm_matrix (np.ndarray): The reduced density matrix of the relevant qubits.
    """
    # Set the quantum state formalism to Density Matrix
    ns.set_qstate_formalism(QFormalism.DM)
    
    # Create the required number of qubits
    qubits = ns.qubits.create_qubits(4*n + 4)

    # Generate the GHZ state for the first group of qubits
    operate(qubits[p['a0']], H)
    for i in range(p['a0'], p['a2']-1):
        operate([qubits[p['a0']], qubits[i+1]], CNOT)

    # Generate the GHZ state for the second group of qubits
    operate(qubits[p['a2']], H)
    for i in range(p['a2'], p['A1']-1):
        operate([qubits[p['a2']], qubits[i+1]], CNOT)

    # Apply noise to the qubits
    for i in range(p['a0'], p['A1']):
        apply_noise(qubits[i], noise, q)

    # Calculate the fidelity of the generated GHZ state
    fidelity_ghz = fidelity(
        [qubits[p['a0']], qubits[p['a1']], qubits[p['b0']], qubits[p['b1']]] + qubits[4:n*2],
        GHZ_party,
        squared=True
    )
    
    # Measure qubits from the first group
    sum_1 = 0
    for i in range(p['a0'], p['a2']):
        if i not in [p['a1'], p['b0'], p['b1']]:
            operate(qubits[i], H)
            meas, prob = ns.qubits.measure(qubits[i])
            if meas == 1:
                sum_1 = sum_1 + 1
                
    # Measure qubits from the second group
    sum_2 = 0
    for i in range(p['a2'], p['A1']):
        if i not in [p['a2'], p['a3'], p['b2']]:
            operate(qubits[i], H)
            meas, prob = ns.qubits.measure(qubits[i])
            if meas == 1:
                sum_2 = sum_2 + 1

    # Apply Z operation to qubit a1 if sum_1 is odd
    if sum_1 % 2 == 1:
        operate(qubits[p['a1']], Z)
                
    # Generate random bits r_a2 and r_a3
    r_a2 = random.randint(0, 1)
    r_a3 = random.randint(0, 1)
    sum_2 = sum_2 + r_a2 + r_a3

    # Apply Z operation to qubit b2 if sum_2 is odd
    if sum_2 % 2 == 1:
        operate(qubits[p['b2']], Z)
                
    # Apply Z operation to qubits a2 and a3 based on random bits
    if r_a2 == 1:
        operate(qubits[p['a2']], Z)
                
    if r_a3 == 1:
        operate(qubits[p['a3']], Z)

    # Calculate the fidelity of the anonymous state with the 3-qubit GHZ state
    fidelity_anonymous = fidelity(
        [qubits[p['a1']], qubits[p['b0']], qubits[p['b1']]],
        GHZ_3,
        squared=True
    )

    # Compute the reduced density matrix of the relevant qubits
    dm_matrix = reduced_dm([
        qubits[p['a1']], qubits[p['b0']], qubits[p['b1']],
        qubits[p['a2']], qubits[p['a3']], qubits[p['b2']]
    ])

    return fidelity_anonymous, fidelity_ghz, dm_matrix
