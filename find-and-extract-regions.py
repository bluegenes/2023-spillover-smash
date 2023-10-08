import pysam
import argparse
from collections import defaultdict

def find_and_extract_genomic_regions(bam_file, fasta_file, output_fasta, output_coords):
    # Identify regions from BAM
    min_starts = defaultdict(lambda: float('inf'))
    max_ends = defaultdict(lambda: 0)

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if not read.is_unmapped:
                reference_name = read.reference_name
                min_starts[reference_name] = min(min_starts[reference_name], read.reference_start)
                max_ends[reference_name] = max(max_ends[reference_name], read.reference_end)

    # Extract regions from FASTA
    with pysam.Fastafile(fasta_file) as fasta, open(output_fasta, 'w') as fasta_out, open(output_coords, 'w') as coords_out:
        for chr, start in min_starts.items():
            end = max_ends[chr]
            sequence = fasta.fetch(chr, start, end)
            fasta_out.write(f">{chr}:{start+1}-{end}\n{sequence}\n")  # +1 to make the start coordinate 1-based
            coords_out.write(f"{chr}\t{start + 1}\t{end}\n")  # Write out the 1-based coordinates

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find genomic regions in BAM and extract from FASTA file.")
    parser.add_argument("--bam", required=True, help="Input BAM file.")
    parser.add_argument("--fasta", required=True, help="Input reference FASTA file.")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA file with extracted segments.")
    parser.add_argument("--output-coords", required=True, help="Output coordinates of the extracted segments.")
    args = parser.parse_args()

    find_and_extract_genomic_regions(args.bam, args.fasta, args.output_fasta, args.output_coords)

