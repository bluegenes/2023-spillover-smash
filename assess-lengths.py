import os
import gzip
import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt

def calculate_sequence_lengths(directory):
    # list to store sequence lengths
    seq_lengths = []

    # loop through all files in the directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            # check if file is .fna.gz
            if file.endswith('.fna.gz'):
                with gzip.open(os.path.join(root, file), 'rt') as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        # append length of each sequence to the list
                        seq_lengths.append(len(record.seq))
    return seq_lengths

def plot_histogram(seq_lengths, output_file):
    # create histogram
    plt.hist(seq_lengths, bins=50, alpha=0.5)
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')
    plt.title('Histogram of Sequence Lengths')
    plt.grid(True)
    # save the figure to a file
    plt.savefig(output_file)

def main():
    # create parser
    parser = argparse.ArgumentParser(description="Calculate sequence lengths and plot histogram.")
    parser.add_argument("-d", "--directory", required=True, help="Directory to search for .fna.gz files")
    parser.add_argument("-o", "--output", default="histogram.png", help="Output file name for the histogram")
    args = parser.parse_args()

    # calculate sequence lengths
    seq_lengths = calculate_sequence_lengths(args.directory)

    # plot histogram
    plot_histogram(seq_lengths, args.output)

if __name__ == "__main__":
    main()

