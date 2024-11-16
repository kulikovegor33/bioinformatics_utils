from typing import List, Union, Dict, Tuple

from scripts.script_dna_rna_tools import is_valid_sequence, transcribe, reverse
from scripts.script_dna_rna_tools import complement, reverse_complement
from scripts.script_fastq_filter import gc_content, average_quality
from scripts.script_fastq_filter import read_fastq, write_fastq


def run_dna_rna_tools(
    *args: Union[str, List[str]]
) -> Union[str, List[str], None]:

    """
    Executes specified procedures (transcription,
    reverse, complementarity) for given
    DNA/RNA sequences.

    Args:
        *args: DNA/RNA sequences and
        the procedure to perform.
               The last argument must be a string
               indicating the procedure.

    Returns:
        The result of the procedure execution:
        - If the procedure is "transcribe", returns
          the transcribed sequence.
        - If the procedure is "reverse",
          returns the reversed sequence.
        - If the procedure is "complement",
          returns the complementary sequence.
        - If the procedure is "reverse_complement",
          returns the reverse complementary sequence.

        If there is one sequence, returns
        a string; if more than one, returns a list of strings.
        If there are no sequences, returns
        None.

    Raises:
        ValueError: If the sequence
        is invalid or the procedure is unknown.
    """

    if not args:
        return None

    procedure = args[-1]
    sequences = args[:-1]

    for seq in sequences:
        if not is_valid_sequence(seq):
            raise ValueError(f"Invalid sequence: {seq}")

    if procedure == "transcribe":
        result = [transcribe(seq) for seq in sequences]
    elif procedure == "reverse":
        result = [reverse(seq) for seq in sequences]
    elif procedure == "complement":
        result = [complement(seq) for seq in sequences]
    elif procedure == "reverse_complement":
        result = [reverse_complement(seq) for seq in sequences]
    else:
        raise ValueError(f"Unknown procedure: {procedure}")

    return result if len(result) > 1 else result[0]


def filter_fastq(input_fastq: str, output_fastq: str,
                 gc_bounds: Tuple[float, float] = (0.0, 100.0),
                 length_bounds: Tuple[int, int] = (0, 2**32),
                 quality_threshold: int = 0) -> None:
    """
    Filters sequences in FASTQ format based on specified criteria.

    Parameters:
    - input_fastq (str): Path to the input FASTQ file.
    - output_fastq (str): Path to the output FASTQ file
where filtered sequences will be saved.
    - gc_bounds (Tuple[float, float]): Minimum and 
maximum GC content percentage for filtering.
    - length_bounds (Tuple[int, int]): Minimum and
maximum sequence length for filtering.
    - quality_threshold (int): Minimum average quality
score for filtering.

    Returns:
    - None: The function writes the filtered sequences
to the output file.
    """

    # Check and create directory for saving filtered data
    output_dir = os.path.dirname(output_fastq)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read data from the FASTQ file
    seqs: Dict[str, Tuple[str, str]] = read_fastq(input_fastq)
    filtered_seqs: Dict[str, Tuple[str, str]] = {}

    for name, (sequence, quality) in seqs.items():
        gc: float = gc_content(sequence)
        length: int = len(sequence)
        avg_quality: int = average_quality(quality)

        # Filtering based on conditions
        if (gc_bounds[0] <= gc <= gc_bounds[1] and
                length_bounds[0] <= length <= length_bounds[1] and
                avg_quality >= quality_threshold):
            filtered_seqs[name] = (sequence, quality)

    # Write filtered data to the output file
    write_fastq(filtered_seqs, output_fastq)

if __name__ == "__main__":
    # Usage of the function
    input_file: str = "example_fastq.fastq"
    output_file: str = "output.fastq"

    filter_fastq(input_file, output_file, gc_bounds=(30.0, 70.0), length_bounds=(50, 300), quality_threshold=20)
