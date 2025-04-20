# Bioinformatics Utils

This project implements object-oriented bioinformatics tools for working with DNA, RNA, and protein sequences, as well as filtering FASTQ files using Biopython.

## Table of Contents
- Installation
- Usage
  - Biological Sequence Classes
    - `BiologicalSequence`
    - `NucleicAcidSequence`
    - `DNASequence`
    - `RNASequence`
    - `AminoAcidSequence`
  - FASTQ Filtering
- Examples

## Installation

Ensure you have Python 3.11 or higher and biopython 1.79 or higher installed.

## Usage

### Biological Sequence Classes

This project follows an OOP approach, where different types of biological sequences are implemented as classes.

#### `BiologicalSequence`
Abstract base class for all biological sequences.
- Supports sequence length (`len(seq)`), indexing (`seq[i]`), and string conversion (`str(seq)`).
- Defines the abstract method `check_alphabet()` for sequence validation.

#### `NucleicAcidSequence`
Intermediate class for nucleic acid sequences (DNA & RNA).
- Implements:
  - `complement()` - Returns the complementary sequence.
  - `reverse()` - Returns the reversed sequence.
  - `reverse_complement()` - Returns the reverse complementary sequence.
- Uses `valid_chars` and `complement_map` to differentiate DNA and RNA.

#### `DNASequence`
Specialized class for DNA sequences, inheriting from `NucleicAcidSequence`.
- Implements `transcribe()` to convert DNA into RNA.

#### `RNASequence`
Specialized class for RNA sequences, inheriting from `NucleicAcidSequence`.

#### `AminoAcidSequence`
Class for protein sequences, inheriting from `BiologicalSequence`.
- Implements `molecular_weight()` to compute protein mass using Biopython.

### FASTQ Filtering

`filter_fastq()` filters sequences from a FASTQ file based on:
- GC content range (`gc_bounds`)
- Sequence length range (`length_bounds`)
- Minimum average quality score (`quality_threshold`)

Uses `SeqIO` and `SeqRecord` from Biopython.

## Examples

### Creating and Using DNA Sequences
```python
dna = DNASequence("ATGC")
print(dna.complement())  # Output: "TACG"
print(dna.transcribe())  # Output: "AUGC"
```

### Filtering FASTQ Files
```python
filter_fastq("input.fastq", "output.fastq", gc_bounds=(30, 70), length_bounds=(50, 300), quality_threshold=20)
```