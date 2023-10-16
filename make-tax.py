import os,sys
import argparse
import csv
import pandas as pd
import screed

def get_gene_accs(fasta):
    accs = []
    with screed.open(fasta) as f:
        for record in f:
            acc = record.name.split(" ")[0]
            accs.append(acc)
    return accs


def build_geneD(fromfile_csv_list, suppressed_records=[]):
    geneD = {}
    for inF in fromfile_csv_list:
        with open(inF) as csvfile:
            r = csv.DictReader(csvfile)
#            for genome_acc in suppressed_records:
#                geneD[genome_acc] = {"dna": [], "protein": []}
            for row in r:
                genome_acc = row['name'].split(' ')[0]
                gene_accs = get_gene_accs(row["genome_filename"])
                prot_accs = []
                if row["protein_filename"]:
                    prot_accs = get_gene_accs(row["protein_filename"])
                geneD[genome_acc] = {"dna": gene_accs, "protein": prot_accs}
    return geneD

def main(args):

    # read in VMR tsv containing virus metadata for exemplar and additional isolate genomes
    vmr = pd.read_csv(args.vmr_tsv, sep='\t')
    assembly_ident_col = 'GenBank Assembly ID'
    genbank_ident_col = 'Virus GENBANK accession'
    refseq_ident_col = 'Virus REFSEQ accession'

    # handle GenBank Failures and suppressed records
    vmr.loc[vmr['GenBank Assembly ID'].isin(args.suppressed_records), 'GenBank Failures'] = 'suppressed'
    vmr['VMR_Accession'] = 'VMR_MSL38_' + vmr['Sort'].astype(str)
    curate_info = ['suppressed', 'no_assembly', 'multiple_acc']#, 'parentheses']#, 'retrieval']

    vmr = vmr.rename(columns={assembly_ident_col: 'ident', genbank_ident_col: 'genbank_ident', refseq_ident_col: 'refseq_ident', 'Virus name(s)': 'name', 'Exemplar or additional isolate': 'exemplar_or_additional'})
    # if GenBank Failures is in curate_info, set ident column to VMR accession
    vmr.loc[vmr['GenBank Failures'].isin(curate_info), 'ident'] = vmr['VMR_Accession']
    vmr['superkingdom'] = 'Viruses'
    vmr = vmr.rename(columns=str.lower)
    lineage_columns = ['superkingdom', 'realm', 'subrealm', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass',
                 'order', 'suborder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'name']
    # some names have ';' in them, e.g. 'invertebrate iridescent virus 6; Chilo iridescent virus'
    # some (one) have ',' in them. Replace that with '_' to prevent issues.
    # Options:
    # - could keep tsv instead of csv, using ',' to sep. Problem: all sourmash taxonomy files are csv, not designed to work with tsv
    # - alternatively, replace ';' with '__' and ',' with '_' in the name column
    vmr['name'] = vmr['name'].str.replace(',', '_')
    vmr['name'] = vmr['name'].str.replace(';', '__')
    vmr['lineage'] = vmr[lineage_columns].fillna('').apply(lambda x: ','.join(x.astype(str)), axis=1)
    vmr = vmr[vmr['ident'].notnull()]
    # build dictionary of genome assembly accession to gene/protein accessions
    geneD = build_geneD(args.fromfile, suppressed_records=args.suppressed_records)
    vmr['gene_accs'] = vmr['ident'].apply(lambda x: ','.join(geneD[x]['dna']))
    vmr['protein_accs'] = vmr['ident'].apply(lambda x: ','.join(geneD[x]['protein']))
    print(vmr.shape)
    tax_columns = ['ident', 'genbank_ident', 'refseq_ident'] + lineage_columns + ['lineage', 'exemplar_or_additional', 'gene_accs', 'protein_accs']
    tax_info = vmr[tax_columns]
    tax_info.to_csv(args.output, sep=',', index=False)

def cmdline(sys_args):
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--fromfile', nargs='+', help='sourmash fromfile csv(s) with name,genome_filename,protein_filename columns', default=["output.vmr/vmr_MSL38_v1.fromfile.csv"])
    p.add_argument('-v', '--vmr-tsv', help='VMR tsv with genbank assembly accessions', default="inputs/VMR_MSL38_v1.acc.tsv")
    p.add_argument('-o', '--output', help='Output taxonomy CSV file', default="output.vmr/vmr_MSL38_v1.taxonomy.csv")
    p.add_argument('-s', '--suppressed-records', nargs='+', help='Suppressed records', default=['GCF_002987915.1', 'GCF_002830945.1', 'GCF_002828705.1', 'GCA_004789135.1'])
    args = p.parse_args()
    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
