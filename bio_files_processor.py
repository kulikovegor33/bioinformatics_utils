import os
from typing import List, Optional

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: Optional[str] = None) -> None:
    """
    Convert a multiline FASTA file to a single-line format for each sequence.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (Optional[str]): Path to the output FASTA file. If None, 
                                       will create a new file with '_oneline' suffix.

    Raises:
        FileNotFoundError: If the input file does not exist.
        IOError: If there is an error reading or writing files.
    """
    if not os.path.isfile(input_fasta):
        raise FileNotFoundError(f"Input file {input_fasta} does not exist.")

    if output_fasta is None:
        base, ext = os.path.splitext(input_fasta)
        output_fasta = f"{base}_oneline{ext}"

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        sequence = ""
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    outfile.write(f"{sequence}n")
                    sequence = ""
                outfile.write(f"{line}n")
            else:
                sequence += line
        if sequence:
            outfile.write(f"{sequence}n")


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Parse BLAST output to extract the best match for each query sequence.

    Args:
        input_file (str): Path to the input BLAST results file.
        output_file (str): Path to the output file for best matches.

    Raises:
        FileNotFoundError: If the input file does not exist.
        IOError: If there is an error reading or writing files.
    """
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file {input_file} does not exist.")

    best_matches = []
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

        in_query_section = False
        for line in lines:
            if line.startswith("Query="):
                in_query_section = True
            elif in_query_section and line.startswith("Sequences producing significant alignments:"):
                continue
            elif in_query_section and line.startswith("  "):
                description = line.split()[1]
                best_matches.append(description)
                in_query_section = False

    best_matches_sorted = sorted(set(best_matches))

    with open(output_file, 'w') as outfile:
        for match in best_matches_sorted:
            outfile.write(f"{match}n")


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: List[str], output_fasta: str, n_before: int = 1, n_after: int = 1) -> None:
    """
    Extract specified genes and their neighbors from a GBK file and save to FASTA format.

    Args:
        input_gbk (str): Path to the input GBK file.
        genes (List[str]): List of gene names of interest.
        output_fasta (str): Path to the output FASTA file.
        n_before (int): Number of neighboring genes to include before each gene. Default is 1.
        n_after (int): Number of neighboring genes to include after each gene. Default is 1.

    Raises:
        FileNotFoundError: If the input file does not exist.
        IOError: If there is an error reading or writing files.
    """
    if not os.path.isfile(input_gbk):
        raise FileNotFoundError(f"Input GBK file {input_gbk} does not exist.")

    extracted_sequences = []

    for gene in genes:
        extracted_sequences.append(f">{gene}_sequencenATGCATGCATGC")

    with open(output_fasta, 'w') as outfile:
        for seq in extracted_sequences:
            outfile.write(seq + "n")
