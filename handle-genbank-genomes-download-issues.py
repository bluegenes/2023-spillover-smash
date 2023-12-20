import os
import argparse

def main(args):
    # Define dirs for info.csv and genomic.fna.gz files
    info_directory = args.info_directory
    genome_directory = args.genome_directory

    # Get the list of files in each directory
    info_files = os.listdir(info_directory)
    genome_files = os.listdir(genome_directory)

    # Get all accessions for which we have a genome file
    genome_accs = {os.path.basename(fn).split('_genomic')[0] for fn in genome_files}

    n_to_delete = 0
    # Iterate through the info file accessions to find those that do not have a corresponding genome file
    for info_file in info_files:
        info_acc = os.path.basename(info_file).split('.info')[0]
        if info_acc not in genome_accs:
            n_to_delete += 1
            # build full path
            info_file_path = os.path.join(info_directory, info_file)
            if args.delete_files:
                os.remove(info_file_path)
                if args.verbose:
                    print(f"Deleted: {info_file_path}")
            elif args.verbose:
                print(f"Would delete: {info_file_path}")
    
    if not args.delete_files:
        print(f"Found {n_to_delete} files to delete.")
    else:
        print(f"Deleted {n_to_delete} files.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find and optionally delete info.csv files without corresponding genomic.fna.gz files.")
    parser.add_argument("--info-directory", default="genbank/info", help="Path to the directory containing info.csv files.")
    parser.add_argument("--genome-directory", default="genbank/genomes", help="Path to the directory containing genomic.fna.gz files.")
    parser.add_argument("--delete-files", action="store_true", help="Actually delete info.csv files without corresponding genome files.")
    parser.add_argument("--verbose", action="store_true", help="Print more information.")
    args = parser.parse_args()
    main(args)





