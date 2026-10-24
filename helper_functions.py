import math

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """

    # Calculating hamming distance and length ignoring indels
    hamming = 0
    length = 0
    for i in range(len(reference_sequence)):
        if reference_sequence[i] == '-' or distant_sequence[i] == '-':
            continue
        elif reference_sequence[i] != distant_sequence[i]:
            hamming += 1

        length += 1

    # Calculate Jukes-Cantor correction
    correction = -3/4*math.log(1-(4/3)*(hamming/length))

    # Return corrected distance
    return (correction*length)


def check_transition(nuc_1, nuc_2):
    if (nuc_1 == 'A' and nuc_2 == 'G') or (nuc_1 == 'G' and nuc_2 == 'A'):
        return True
    elif (nuc_1 == 'C' and nuc_2 == 'T') or (nuc_1 == 'T' and nuc_2 == 'C'):
        return True

    return False

def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """

    # K = -1/2 ln[1 - 2P - Q)] - 1/4 ln[1 - 2Q]
    # P is the fraction of nucleotides that a simple count reveals to be transitions, and
    # Q is the fraction of nucleotides that a simple count reveals to be transversions

    # Calculting P and Q counts
    p_count = 0
    q_count = 0
    length = 0

    for i in range(len(reference_sequence)):
            if reference_sequence[i] == '-' or distant_sequence[i] == '-':
                continue
            elif reference_sequence[i] != distant_sequence[i]:
                if check_transition(reference_sequence[i], distant_sequence[i]):
                    p_count += 1
                else:
                    q_count += 1
    
            length += 1

    p = p_count/length
    q = q_count/length
    
    correction = -1/2*math.log(1-2*p-q) - 1/4*math.log(1-2*q)

    return (correction*length)
