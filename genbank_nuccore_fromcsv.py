import sys
import gzip
import os
import csv
import argparse
from Bio import Entrez, SeqIO
from Bio.Seq import UndefinedSequenceError

# Set your email address (required by NCBI)
Entrez.email = "ntpierce@ucdavis.edu"

def download_and_join_fastas(accession_list, nucleotide_file, protein_file):
    """
    Downloads nucleotide and/or protein sequences for the given GenBank accession number
    and saves them as gzip-compressed FASTA files with the specified filenames, or with
    default filenames based on the accession number.
    
    Args:
        accession (str): The GenBank accession number to download.
        nucleotide_file (str): The filename for the nucleotide FASTA file.
        protein_file (str): The filename for the protein FASTA file.

    Returns:
        tuple: A tuple containing (accession, organism_name, nucleotide_file, protein_file).
    """
    # Use Entrez to download the record
    nucl_seqs = []
    prot_seqs = []
    failed_accinfo = []
    for acc in accession_list:
        acname = ""
        if ':' in acc: # e.g. RNA1:KT601119
            acname, acc = acc.split(':')
            acname = acname + ':'
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            # organism_name = record.annotations["organism"]

            # Nucleotide sequence
            nucleotide_sequence = None
            try:
                nucleotide_sequence = str(record.seq)
                nucl_seqs.append((acname + acc, nucleotide_sequence))
            except UndefinedSequenceError:
                sys.stderr.write(f"Could not find sequence for accession {acc}: record.seq is undefined.\n")
                failed_accinfo.append((acname + acc, "DNA"))

            # Protein sequences
            protein_sequences = []
            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    protein_sequences.append((feature.qualifiers["translation"][0], feature.qualifiers["protein_id"][0]))
            if not protein_sequences:
                failed_accinfo.append((acname + acc, "protein"))
            prot_seqs.append((acname + acc, protein_sequences))

        except Exception as e:
            sys.stderr.write(f"Failed to download FASTA for {acc}: {str(e)}\n")
            failed_accinfo.append((acname + acc, "DNA"))
            failed_accinfo.append((acname + acc, "protein"))
            
    # Save the sequences to gzip-compressed FASTA files
    if nucl_seqs:
        with gzip.open(nucleotide_file, "wt") as f:
            for ac, seq in enumerate(nucl_seqs):
                f.write(f">{ac}\n{seq}\n")
    if prot_seqs:
        with gzip.open(protein_file, "wt") as f:
            for ac, seqs in prot_seqs:
                for i, (protein_sequence, protein_id) in enumerate(seqs):
                    f.write(f">{ac}_CDS{i+1}|{protein_id}\n{protein_sequence}\n")
    else:
        protein_file = None

    return failed_accinfo, nucleotide_file, protein_file if prot_seqs else None

def main(args):
    """
    Processes the CSV file and downloads nucleotide sequences for the given GenBank accession numbers.
    Also writes a single fileinfo CSV file with details for all accessions.
    
    Args:
        args (Namespace): The command-line arguments.

    Returns:
        None
    """
    all_fileinfo = []
    failed_accs = []
    outdir = args.output_fastadir

    with open(args.csvfile, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row['name']
            ident = name.split(' ')[0]
            accessions = row['genbank_accessions'].split(';')
            nucleotide_file = f"{outdir}/{ident}.fna.gz"
            protein_file = f"{outdir}/{ident}.faa.gz"
            failed_accinfo, nucl_file, prot_file = download_and_join_fastas(accessions, nucleotide_file, protein_file)
            failed_accs.append((name, *failed_accinfo))
            all_fileinfo.append((name, nucl_file, prot_file))
    
    # Write fileinfo to a CSV file
    with open(args.fileinfo, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["ident", "name", "genome_file", "protein_file"])
        for info in all_fileinfo:
            writer.writerow(info)
    
    # write failed accessions to a CSV file
    with open(args.failed, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["name", "accession", "sequence_type"])
        for info in failed_accs:
            writer.writerow(info)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(description="Download sequences for a list of GenBank accessions from a CSV file")
    p.add_argument("csvfile", type=str, help="The CSV file containing accession numbers and download filenames")
    p.add_argument("--fileinfo", type=str, help="The filename for the fileinfo CSV that contains details for all accessions")
    p.add_argument('-o','--output-fastadir', type=str, help="The directory to save the output FASTA files", default=".")
    p.add_argument("--failed", type=str, help="The filename for the failed accessions CSV file")

    # Parse command line arguments
    args = p.parse_args()

    # Process the CSV file and download sequences
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
