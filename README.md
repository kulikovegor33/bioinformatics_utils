# DNA/RNA Tools and FASTQ Filter

This project contains functions for working with DNA/RNA sequences and filtering sequences in FASTQ format.

## Table of Contents

- Installation
- Usage
  - Function run_dna_rna_tools
  - Function filter_fastq
  - Function convert_multiline_fasta_to_oneline
  - Function parse_blast_output
  - Function select_genes_from_gbk_to_fasta
- Exceptions
- Examples

## Installation

Ensure you have Python 3.6 or higher installed. Then, install the required libraries if they are not already installed.

```bash
pip install -r requirements.txt
```
## Usage

### Function run_dna_rna_tools

This function performs specified procedures (transcription, reverse, complement) for given DNA/RNA sequences.

```python
def run_dna_rna_tools(*args: Union[str, List[str]]) -> Union[str, List[str], None]:
```
#### Arguments

- *args: DNA/RNA sequences and the procedure to perform. The last argument should be a string indicating the procedure.

#### Return Value

- Result of the procedure:
  - If the procedure is "transcribe", returns the transcribed sequence.
  - If the procedure is "reverse", returns the reversed sequence.
  - If the procedure is "complement", returns the complementary sequence.
  - If the procedure is "reverse_complement", returns the reverse complementary sequence.

If there is one sequence, it returns a string; if there are multiple, it returns a list of strings. If there are no sequences, it returns None.

### Function filter_fastq

```python
filter_fastq(input_fastq: str, output_fastq: str, 
              gc_bounds: Tuple[float, float] = (0.0, 100.0), 
              length_bounds: Tuple[int, int] = (0, 2**32), 
              quality_threshold: int = 0) -> None
```
#### Parameters:

- input_fastq (str): Path to the input FASTQ file.
- output_fastq (str): Path to the output FASTQ file where filtered sequences will be saved.
- gc_bounds (Tuple[float, float]): Minimum and maximum GC content for filtering.
- length_bounds (Tuple[int, int]): Minimum and maximum sequence length for filtering.
- quality_threshold (int): Minimum average quality score for filtering.

### Function convert_multiline_fasta_to_oneline

This function converts multi-line FASTA files into single-line format.

```python
def convert_multiline_fasta_to_one_line(fasta_file: str) -> Dict[str, str]:
```
#### Arguments

- fasta_file: Path to the multi-line FASTA file.

#### Return Value

- A dictionary with sequence identifiers as keys and single-line sequences as values.

### Function parse_blast_output

This function parses BLAST output files and extracts relevant information.

```python
def parse_blast_output(blast_file: str) -> List[Dict[str, Union[str, float]]]:
```
#### Arguments

- blast_file: Path to the BLAST output file.

#### Return Value

- A list of dictionaries containing parsed BLAST results.

### Function select_genes_from_gbk_to_fasta

This function selects specific genes from GenBank files and converts them into FASTA format.

```python
def select_genes_from_gbk_to_fasta(gbk_file: str, gene_list: List[str]) -> Dict[str, str]:
```
#### Arguments

- gbk_file: Path to the GenBank file.
- gene_list: List of gene names to extract.

#### Return Value

- A dictionary with gene names as keys and their corresponding FASTA sequences as values.

## Exceptions

- ValueError: If a sequence is invalid or an unknown procedure is specified.

## Examples

### Example usage of function run_dna_rna_tools

```python
result = run_dna_rna_tools("ACGT", "transcribe")
print(result)  # Output: "UGCA"
```

### Example usage of function filter_fastq

```python
if name == "main":
    input_file = "example_fastq.fastq"
    output_file = "output.fastq"

    filter_fastq(input_file, output_file, gc_bounds=(30.0, 70.0), length_bounds=(50, 300), quality_threshold=20)
```

### Example usage of function convert_multiline_fasta_to_oneline

```python
sequences = convert_multiline_fasta_to_oneline("input.fasta")
print(sequences)  # Output: {'seq1': 'ATCG...', 'seq2': 'GCTA...'}
```
### Example usage of function parse_blast_output

```python
results = parse_blast_output("output.blast")
print(results)  # Output: [{'query': 'seq1', 'subject': 'subject1', 'score': 100.0}, ...]
```

### Example usage of function select_genes_from_gbk_to_fasta

```python
fasta_sequences = select_genes_from_gbk_to_fasta("genome.gbk", ["gene1", "gene2"])
print(fasta_sequences)  # Output: {'gene1': 'ATCG...', 'gene2': 'GCTA...'}
```
