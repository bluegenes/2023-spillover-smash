import sys
import gzip
import csv
import argparse


def read_input_fileinfo(fileinfo):
    """Reads in fileinfo csv file and returns a list of tuples containing
    the VMR accession number, the path to the nucleotide fasta file, and the
    path to the protein fasta file (optional)."""
    nucl, prot = [],[]
    for fileinfo_csv in fileinfo:
        with open(str(fileinfo_csv)) as inF:
            reader = csv.reader(inF)
            for row in reader:
                # row = name,genome,protein
                nucl.append(row[1])
                if row[2]:
                    prot.append(row[2])
    sys.stderr.write(f"nucl: {nucl}\n")
    sys.stderr.write(f"prot: {prot}\n")
    return nucl, prot


def combine_fasta_files(fastas, outF):
    """Combines fasta files into a single file."""
    with gzip.open(str(outF), 'wt') as outF:
        for fasta in fastas:
            with gzip.open(str(fasta), 'rt') as inF:
                for line in inF:
                    outF.write(line)


def main(args):
    nucl, prot= read_input_fileinfo(args.input_fileinfo)
    # combine fasta files
    combine_fasta_files(nucl, args.nucl_out)
    if prot:
        combine_fasta_files(prot, args.prot_out)

    # write fileinfo
    with open(str(args.fileinfo_out), "w") as outF:
        if prot:
            outF.write(f"{args.curated_acc},{args.nucl_out},{args.prot_out}\n")
        else:
            outF.write(f'{args.curated_acc},{args.nucl_out},\n')


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(description="Combines fasta files from multiple sources into a single file.")
    p.add_argument("--input-fileinfo", nargs="+",
        help="File containing information about fasta files to combine. " +
        "Each line should contain the VMR accession number, the path to the " +
        "nucleotide fasta file, and the path to the protein fasta file " +
        "(optional).")
    p.add_argument("--curated-acc", type=str, help="VMR accession number for the combined fasta file.")
    p.add_argument("--nucl-out", help="Path to output nucleotide fasta file.")
    p.add_argument("--prot-out", help="Path to output protein fasta file.")
    p.add_argument("--fileinfo-out", help="Path to output fileinfo file.")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)