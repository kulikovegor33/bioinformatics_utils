import os
from typing import Dict, Tuple

def gc_content(sequence: str) -> float:
    """Calculates the guanine-cytosine (GC) content of a sequence.

    Args:
        sequence (str): DNA sequence.

    Returns:
        float: Percentage of GC content in the sequence.
    """
    g_count = sequence.count('G') + sequence.count('g')
    c_count = sequence.count('C') + sequence.count('c')
    total_count = len(sequence)

    if total_count == 0:
        return 0.0

    return (g_count + c_count) / total_count * 100

def average_quality(quality: str) -> float:
    """Calculates the average quality of a sequence based on ASCII codes.

    Args:
        quality (str): Quality string (usually in Phred format).

    Returns:
        float: Average quality value.
    """
    return sum(ord(q) - 33 for q in quality) / len(quality)

def read_fastq(input_fastq: str) -> Dict[str, Tuple[str, str]]:
    """Reads a FASTQ file and returns its contents.

    Args:
        input_fastq (str): Path to the input FASTQ file.

    Returns:
        Dict[str, Tuple[str, str]]: Dictionary with sequence names as keys
        and tuples (sequence, quality) as values.
    """
    seqs = {}
    with open(input_fastq, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            f.readline()  # Skip the '+' line
            quality = f.readline().strip()
            seqs[header] = (sequence, quality)
    return seqs

def write_fastq(seqs: Dict[str, Tuple[str, str]], output_fastq: str):
    """Writes filtered sequences to a FASTQ file.

    Args:
        seqs (Dict[str, Tuple[str, str]]): Filtered sequences.
        output_fastq (str): Path to the output FASTQ file.
    Raises:
        FileExistsError: If the output file already exists.
    """
    if os.path.exists(output_fastq):
        raise FileExistsError(f"File {output_fastq} already exists. Please choose a different name.")
    with open(output_fastq, 'w') as f:
        for name, (sequence, quality) in seqs.items():
            f.write(f"{name}n{sequence}n+n{quality}n")
