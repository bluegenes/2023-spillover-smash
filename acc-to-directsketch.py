import argparse
import csv


def extract_accs(id_col):
    "split muliple NCBI accs; split names from segment accessions"
    split_ids = []
    if id_col:
        ids = id_col.split(";")
        for id in ids:
            if '(' in id:
                acc_only = id.split('(')[0]
                acc_only = acc_only.strip()
            if ':' in id:
                acc_only = id.split(':')[1]
                acc_only = acc_only.strip()
                split_ids.append(acc_only)
            else:
                id = id.strip()
                split_ids.append(id)
    return split_ids

def main(args):
    acc2info = {}
    bad_accs = set()
    
    # Read good accessions
    with open(args.good_acc, 'rt') as gacc:
        r = csv.reader(gacc, delimiter='\t')
        for row in r:
            if row[0].startswith('#'):
                continue
            acc = row[0]
            gcf_acc = acc.replace("GCA", "GCF")
            genome_length = row[25]
            ftp_path = row[19]
            acc2info[acc] = (genome_length, ftp_path)
            acc2info[gcf_acc] = (genome_length, ftp_path)
    
    # Read bad accessions
    with open(args.bad_acc, 'rt') as bacc:
        r = csv.reader(bacc, delimiter='\t')
        for row in r:
            if row[0].startswith('#'):
                continue
            acc = row[0]
            gcf_acc = acc.replace("GCA", "GCF")
            bad_accs.add(acc)
            bad_accs.add(gcf_acc)
    
    # Process VMR file
    with open(args.vmr_file, 'rt') as infp, \
         open(args.ds_csv, 'w') as ds_csv, \
         open(args.curated_ds, 'w') as curated_ds, \
         open(args.curate_info, 'w') as curate_info, \
         open(args.suppressed, 'w') as suppressed, \
         open(args.lengths, 'w') as lengths:
        
        ds_csv.write("accession,name,ftp_path\n")
        curated_ds.write("accession,name,moltype,md5sum,download_filename,url,range\n")
        curate_info.write("name,assembly_failure,genbank_accessions\n")
        lengths.write("accession,length\n")
        
        r = csv.DictReader(infp, delimiter='\t')
        suppressed.write(','.join(r.fieldnames) + '\n')
        
        for row in r:
            acc = row['GenBank Assembly ID']
            virus_name = f"{row['Virus name(s)']} {row['Virus isolate designation']}".strip().replace(',', ';')
            
            if acc and acc not in bad_accs:
                name = f"{acc} {virus_name}"
                genome_length, ftp_path = acc2info.get(acc, ("", ""))
                ds_csv.write(f"{acc},{name},{ftp_path}\n")
                lengths.write(f"{acc},{genome_length}\n")
            else:
                if acc:
                    print(f"Skipping {acc} for {virus_name} due to historical suppression or failure.")
                    suppressed.write(','.join(row.values()) + '\n')

                gb_col = row["Virus GENBANK accession"]
                if gb_col == "":
                    continue
                vmr_acc = f"{args.basename}_{row['Species Sort']}"
                curated_name = f"{vmr_acc} {virus_name}".strip()
                dl_filename = f"genbank/curated/{vmr_acc}.fna"
                # here, handle multiple acc input styles:
                gb_acc = extract_accs(gb_col)
                
                dl_links = []
                for gba in gb_acc:
                    if gba:
                        dl_links.append(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={gba}&rettype=fasta&retmode=text")
                dl_info = ";".join(dl_links)
                
                range = ""
                # if there's a range provided, add it.
                # NOTE: accs with range should not have multiple accs to join
                if "(" in gb_col and ")" in gb_col:
                    range = gb_col[gb_col.find("("):gb_col.find(")")+1]
                    range = range.replace('.', ':')
                # write directsketch download file
                if dl_info is not None:
                    curated_ds.write(f"{vmr_acc},{curated_name},DNA,,{dl_filename},{dl_info},{range}\n")
                else:
                    print("dl info was None")
                # write info on failure reason
                fail_reason = row['GenBank Failures']
                curate_info.write(f"{curated_name},{fail_reason},{gb_acc}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GenBank accession data and generate directsketch files.")
    parser.add_argument("--vmr_file", required=True, help="Input VMR file.")
    parser.add_argument("--good_acc", required=True, help="Path to the current GenBank assembly summary file.")
    parser.add_argument("--bad_acc", required=True, help="Path to the historical (suppressed) GenBank assembly summary file.")
    parser.add_argument("--ds_csv", required=True, help="Output directsketch CSV file.")
    parser.add_argument("--curated_ds", required=True, help="Output curated directsketch CSV file.")
    parser.add_argument("--curate_info", required=True, help="Output curation information CSV file.")
    parser.add_argument("--suppressed", required=True, help="Output suppressed accessions CSV file.")
    parser.add_argument("--lengths", required=True, help="Output genome lengths CSV file.")
    parser.add_argument("--basename", required=True, help="Base name for output records.")
    args = parser.parse_args()
    main(args)

