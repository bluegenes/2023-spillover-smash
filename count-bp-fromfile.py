import csv
import screed
import argparse
import os

def count_sequences(filename):
    """ Count the number of basepairs or amino acids in a file. """
    count = 0
    if os.path.exists(filename):
        with screed.open(filename) as file:
            for record in file:
                count += len(record.sequence)
    else:
        count = None
    return count

def process_csv(input_csv, output_csv):
    """ Process the CSV file to count basepairs and amino acids, and write to a new CSV. """
    with open(input_csv, 'r') as infile:
        total_rows = sum(1 for _ in csv.reader(infile)) - 1  # Subtract 1 for header

    progress_interval = total_rows // 20  # 5% of total rows
    if progress_interval == 0:
        progress_interval = 1

    print(f"starting length counting for {total_rows} rows in '{input_csv}'")
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['basepair_count', 'amino_acid_count']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        writer.writeheader()
        for i, row in enumerate(reader, start=1):
            genome_file = row['genome_filename']
            protein_file = row.get('protein_filename', '')  # Handle missing protein_filename
            basepair_count = count_sequences(genome_file)
            amino_acid_count = count_sequences(protein_file) if protein_file else None
            row['basepair_count'] = basepair_count
            row['amino_acid_count'] = amino_acid_count
            writer.writerow(row)

            if i % progress_interval == 0 or i == total_rows:
                print(f"Processed {i}/{total_rows} rows ({i/total_rows*100:.2f}%)")

def main():
    parser = argparse.ArgumentParser(description="Count basepairs and amino acids from CSV files.")
    parser.add_argument("input_csv", help="Input CSV file with genome and protein filenames")
    parser.add_argument("output_csv", help="Output CSV file to write the counts")
    args = parser.parse_args()

    process_csv(args.input_csv, args.output_csv)

if __name__ == "__main__":
    main()

