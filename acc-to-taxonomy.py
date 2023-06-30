#! /usr/bin/env python
import sys
import gzip
import os
import csv
import argparse

from Bio import Entrez, SeqIO

def get_taxonomic_lineage(accession):
    # Fetch the record for the given accession in XML format
    handle = Entrez.efetch(db="nucleotide", id=accession, retmode="xml")
    record = Entrez.read(handle)
    handle.close()

    # Extract the taxonomic lineage from the record
    taxonomic_lineage = None
    for feature in record[0]["GBSeq_feature-table"]:
        if feature["GBFeature_key"] == "source":
            qualifiers = feature["GBFeature_quals"]
            for qualifier in qualifiers:
                if qualifier["GBQualifier_name"] == "db_xref" and qualifier["GBQualifier_value"].startswith("taxon:"):
                    taxonomic_lineage = qualifier["GBQualifier_value"].split(":")[1]
                    break

    # Fetch the taxonomic information using the taxon ID
    if taxonomic_lineage:
        handle = Entrez.efetch(db="taxonomy", id=taxonomic_lineage, retmode="xml")
        taxonomy_record = Entrez.read(handle)
        handle.close()

        # Extract the taxonomic lineage with ranks
        lineage = taxonomy_record[0]["LineageEx"]
        lineageD = {}
        for taxon in lineage:
            rank = taxon.get("Rank")
            name = taxon.get("ScientificName")
            if rank and name:
                lineageD[rank] = name
        lineageD["name"] = taxonomy_record[0]["ScientificName"] 
    return lineageD

def main(args):
    # Set your email address (required by NCBI)
    Entrez.email = args.email
    # Set the ICTV ranks to retrieve
    ncbi_ranks = ['superkingdom', "clade", "subrealm", "kingdom", "subkingdom", "phylum", "subphylum", \
                  "class", "subclass", "order", "suborder", "family", "subfamily", \
                    "genus", "subgenus", "species", "subspecies", "name"]

    # load accessions from file if specified, otherwise use the single accession
    if args.from_csv:
        with open(args.from_csv, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            accessions = [row['AccessionNumber'].lstrip() for row in reader]
    elif args.from_file:
        with open(args.from_file, "r") as acc_file:
            accessions = acc_file.read().split("\n")
            print(accessions)
    elif args.accession:
        accessions = [args.accession]
    else:
        sys.exit("No accession specified")
    # open output file to write to:
    with open(args.output, "w") as csvfile:
        fields = ["accession"] + ncbi_ranks
        w = csv.DictWriter(csvfile, fieldnames=fields)
        # write header row
        w.writeheader()
        for acc in accessions:
            taxonomic_lineage = {}
            # Fetch the record for the given accession
            taxonomic_lineage = get_taxonomic_lineage(acc)
            taxonomic_lineage['accession'] = acc
            w.writerow(taxonomic_lineage)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    parser = argparse.ArgumentParser(description="Download nucleotide and/or protein sequences for a GenBank accession")
    parser.add_argument('--accession', type=str, help="The GenBank accession for which to retrieve tax information")
    parser.add_argument('--from-file', type=str, help="The filename for the list of accessions to retrieve tax information for")
    parser.add_argument('--from-csv', type=str, help="The filename for the csv including accessions to retrieve tax information for")
    parser.add_argument('--email', type=str, help="Your email address (required by NCBI)", default="ntpierce@ucdavis.edu")
    parser.add_argument('-o', '--output', type=str, help="The filename for the output CSV file", default="spillover.genbank-lineages.csv")
    # Parse command line arguments
    args = parser.parse_args()

    # Download sequences for the specified accession
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)