import sys
import argparse
import pandas as pd
from Bio import Entrez

# Set email address for NCBI access
Entrez.email = "ntpierce@gmail.com"

def retrieve_assembly_accession(identifier):
    """
    Retrieve the NCBI Assemblies dataset identifier for the given RefSeq or GenBank identifier.

    Args:
        identifier (str): The RefSeq or GenBank identifier.

    Returns:
        str: The NCBI Assemblies dataset identifier, or None if no Assembly is found for the given identifier.
    """
    handle = Entrez.esearch(db="assembly", term=identifier)
    record = Entrez.read(handle)
    handle.close()
    if record["Count"] == "0":
        return None
    else:
        assembly_id = record["IdList"][0]
        handle = Entrez.efetch(db="assembly", id=assembly_id, rettype="docsum")
        record = Entrez.read(handle)
        handle.close()
        dataset_id = record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
        return dataset_id

def extract_accs(ids):
    "split names from segment accessions"
    split_ids = []
    for id in ids:
        if ':' in id:
            acc_only = id.split(':')[1]
            acc_only = acc_only.strip()
            split_ids.append(acc_only)
        else:
            return ids
    return split_ids

def find_assembly_accessions(row):
    """
    Adds two columns to the input pandas DataFrame containing the RefSeq and GenBank Assembly IDs
    corresponding to the RefSeq and GenBank IDs in the given row.

    Parameters:
    -----------
    row : pandas Series
        A pandas Series containing the RefSeq and GenBank IDs for a single row
        of a pandas DataFrame. The RefSeq ID should be in the "RefSeq ID" column
        and the GenBank ID should be in the "GenBank ID" column.

    Returns:
    --------
    None
    """
    # Get the RefSeq and GenBank IDs from the row
    refseq_ids, genbank_ids = [],[]
    refseq_assembly_id, genbank_assembly_id = "",""
    print(row["Virus REFSEQ accession"])
    if not pd.isnull(row["Virus REFSEQ accession"]):
        refseq_ids = row["Virus REFSEQ accession"].split(";")
        refseq_ids = extract_accs(refseq_ids)
        print(refseq_ids)
    print(row["Virus GENBANK accession"])
    if not pd.isnull(row["Virus GENBANK accession"]):
        genbank_ids = row["Virus GENBANK accession"].split(";")
        genbank_ids = extract_accs(genbank_ids)
    # Use the first ID in the list to look up the assembly ID
    if refseq_ids and "(" not in refseq_ids[0] and ")" not in refseq_ids[0]:
        refseq_assembly_id = retrieve_assembly_accession(refseq_ids[0])
    if genbank_ids and "(" not in genbank_ids[0] and ")" not in genbank_ids[0]:
        genbank_assembly_id = retrieve_assembly_accession(genbank_ids[0])
    # If there are multiple IDs in the list, check that they all correspond to the same assembly ID
    if len(refseq_ids) > 1:
        for refseq_id in refseq_ids[1:]:
            if "(" not in refseq_id and ")" not in refseq_id:
                if retrieve_assembly_accession(refseq_id) != refseq_assembly_id:
                    raise ValueError("Multiple RefSeq accessions correspond to different assembly IDs.")
    if len(genbank_ids) > 1:
        for genbank_id in genbank_ids[1:]:
            if "(" not in genbank_id and ")" not in genbank_id:
                if retrieve_assembly_accession(genbank_id) != genbank_assembly_id:
                    raise ValueError("Multiple GenBank accessions correspond to different assembly IDs.")
    # Add the assembly IDs to the row
    row["RefSeq Assembly ID"] = refseq_assembly_id
    row["GenBank Assembly ID"] = genbank_assembly_id
    return row


def main(args):
    # VMR Reference file, v37
    vmr_file = "inputs/VMR_21-221122_MSL37.xlsx"
    vmr = pd.read_excel(vmr_file, sheet_name="VMRb37")
    vInfo = vmr.apply(find_assembly_accessions, axis=1)
    vInfo.to_csv('outputs/VMR_21-221122_MSL37.csv.gz')

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("-x", "--xlsx", default= "inputs/VMR_21-221122_MSL37.xlsx")
    p.add_argument("-o", "--output-csv", default='outputs/VMR_21-221122_MSL37.csv.gz')
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)


def test_retrieve_assembly_accession():
    # Test valid input
    assert retrieve_assembly_accession("MK064563") == "GCF_003950175.1"
    # Test valid multiple inputs:
#    assert get_assembly_id("dsRNA 1: HG975302; dsRNA 2: HG975303; dsRNA 3: HG975304; dsRNA 4: HG975305") == "GCA_013138185.1"
    assert retrieve_assembly_accession("HG975302") == "GCA_013138185.1"
    # Test invalid input
    assert retrieve_assembly_accession("invalid_identifier") == None
