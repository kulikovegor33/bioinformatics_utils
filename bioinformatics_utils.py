import logging
from typing import Tuple, List
from abc import ABC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

class BiologicalSequence(ABC):
    def __init__(self, sequence: str) -> None:
        """Initialize a biological sequence and validate its alphabet."""
        self.sequence: str = sequence.upper()
        if not self.check_alphabet():
            raise ValueError("Invalid sequence alphabet.")

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)
    
    def __getitem__(self, index: int) -> str:
        """Return the character at the index."""
        return self.sequence[index]
    
    def __str__(self) -> str:
        """Return the string representation of the sequence."""
        return self.sequence
    
    def __repr__(self) -> str:
        """Return a representation of the sequence object."""
        return f"{self.__class__.__name__}('{self.sequence}')"

    def check_alphabet(self) -> bool:
        """Check if the sequence contains only valid characters."""
        pass

class NucleicAcidSequence(BiologicalSequence):
    complement_map: dict = {}
    valid_chars: set = set()
    
    def complement(self) -> "NucleicAcidSequence":
        """Return the complementary sequence."""
        return self.__class__(self.sequence.translate(self.complement_map))
    
    def reverse(self) -> "NucleicAcidSequence":
        """Return the reversed sequence."""
        return self.__class__(self.sequence[::-1])
    
    def reverse_complement(self) -> "NucleicAcidSequence":
        """Return the reverse complementary sequence."""
        return self.reverse().complement()
    
    def check_alphabet(self) -> bool:
        """Validate if the sequence consists of allowed nucleotide characters."""
        return set(self.sequence).issubset(self.valid_chars)

class DNASequence(NucleicAcidSequence):
    complement_map: dict = str.maketrans("ATGC", "TACG")
    valid_chars: set = {"A", "T", "G", "C"}
    
    def transcribe(self) -> "RNASequence":
        """Return the transcribed RNA sequence."""
        return RNASequence(self.sequence.replace("T", "U"))

class RNASequence(NucleicAcidSequence):
    complement_map: dict = str.maketrans("AUGC", "UACG")
    valid_chars: set = {"A", "U", "G", "C"}

class AminoAcidSequence(BiologicalSequence):
    valid_chars: set = set("ACDEFGHIKLMNPQRSTVWY")
    
    def check_alphabet(self) -> bool:
        """Validate if the sequence consists of allowed aminoacid characters."""
        return set(self.sequence).issubset(self.valid_chars)
    
    def molecular_weight(self) -> float:
        """Calculate the molecular weight of the aminoacid sequence."""
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        return ProteinAnalysis(self.sequence).molecular_weight()


def filter_fastq(input_fastq: str, output_fastq: str,
                 gc_bounds: Tuple[float, float] = (0, 100),
                 length_bounds: Tuple[int, float] = (0, float('inf')),
                 quality_threshold: int = 0,
                 log_path: str = "filter.log") -> None:
    """Filter sequences in a FASTQ file based on GC content, length, and quality."""
    filtered_records: List[SeqRecord] = []
    total = 0
    passed = 0

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s'
    )

    try:
        for record in SeqIO.parse(input_fastq, "fastq"):
            total += 1
            gc_content: float = gc_fraction(record.seq) * 100
            avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record)

            if (gc_bounds[0] <= gc_content <= gc_bounds[1] and
                length_bounds[0] <= len(record) <= length_bounds[1] and
                avg_quality >= quality_threshold):
                filtered_records.append(record)
                passed += 1
        
        SeqIO.write(filtered_records, output_fastq, "fastq")
        logging.info(f"{passed}/{total} records passed filtering and go to {output_fastq}")
    except Exception as e:
        logging.error(f"Error file {input_fastq}: {e}")
        raise

    