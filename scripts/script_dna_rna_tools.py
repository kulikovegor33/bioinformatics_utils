valid_dna = set("ATGCatgc")
valid_rna = set("AUGCaugc")
complement_map_dna = str.maketrans("ATGCatgc", "TACGtacg")
complement_map_rna = str.maketrans("AUGCaugc", "UACGuacg")


def is_valid_sequence(sequence: str) -> bool:
    """
    Checks if the sequence is a valid DNA or RNA sequence.

    Args:
        sequence (str): The sequence to check.

    Returns:
        bool: True if the sequence is valid, otherwise False.
    """
    return all(base in valid_dna or base in valid_rna for base in sequence)


def transcribe(dna_sequence: str) -> str:
    """
    Transcribes DNA to RNA by replacing 'T' with 'U'.

    Args:
        dna_sequence (str): The DNA sequence.

    Returns:
        str: The transcribed RNA sequence.
    """
    return dna_sequence.replace("T", "U").replace("t", "u")


def reverse(sequence: str) -> str:
    """
    Returns the reverse of the sequence.

    Args:
        sequence (str): The original sequence.

    Returns:
        str: The reversed sequence.
    """
    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Finds the complementary sequence.

    Args:
        sequence (str): The original sequence.

    Returns:
        str: The complementary sequence.
    """
    if "U" in sequence or "u" in sequence:
        return sequence.translate(complement_map_rna)
    else:
        return sequence.translate(complement_map_dna)


def reverse_complement(sequence: str) -> str:
    """
    Finds the reverse complementary sequence.

    Args:
        sequence (str): The original sequence.

    Returns:
        str: The reverse complementary sequence.
    """
    return reverse(complement(sequence))
