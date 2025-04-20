import argparse
from typing import Tuple
from bioinformatics_utils import filter_fastq

def parse_args():
    parser = argparse.ArgumentParser(description="Filter FASTQ file by GC content, length, and quality.")
    parser.add_argument("input_fastq", help="Path to input FASTQ file")
    parser.add_argument("output_fastq", help="Path to output FASTQ file")
    parser.add_argument("--gc_min", type=float, default=0.0, help="Minimum GC content")
    parser.add_argument("--gc_max", type=float, default=100.0, help="Maximum GC content")
    parser.add_argument("--length_min", type=int, default=0, help="Minimum sequence length")
    parser.add_argument("--length_max", type=int, default=float("inf"), help="Maximum sequence length")
    parser.add_argument("--quality", type=int, default=0, help="Minimum average quality score")
    return parser.parse_args()

def main():
    args = parse_args()
    gc_bounds = (args.gc_min, args.gc_max)
    length_bounds = (args.length_min, args.length_max)
    filter_fastq(args.input_fastq, args.output_fastq, gc_bounds, length_bounds, args.quality)

if __name__ == "__main__":
    main()
