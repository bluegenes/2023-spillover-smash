import os, sys
import argparse
import pandas as pd
import numpy as np

def main(args):

    # build output file names
    out_blastn_missed = args.output_base + '.blastn.missed.tsv'
    out_merged_missed = args.output_base + '.merged.missed.tsv'
    out_besthit_classif = args.output_base + '.bestblast.tsv'
    out_besthit_clean = args.output_base + '.bestblast.clean.tsv'

    sDF = pd.read_csv(args.spillover_csv, index_col=0)
    blastn_df = pd.read_csv(args.blastn, sep='\t')
    blastx_df = pd.read_csv(args.blastx, sep='\t')
    existing_spillover_columns = ['IndividualID', 'AccessionNumber', 'Virus', 'VirusGenus', 'VirusSpecies', 'VirusFamily', 'HostSpecies', 'HostGenus', 'HostFamily']

    merged_bn = sDF.merge(blastn_df, on=existing_spillover_columns, how='left')
    merged_bx = sDF.merge(blastx_df, on=existing_spillover_columns, how='left')

    blastn_missed = merged_bn[merged_bn['lineage'].isnull()]
    blastn_missed.to_csv(out_blastn_missed, sep='\t', index=False)

    blastn_classif = merged_bn[merged_bn['lineage'].notnull()]
    print(blastn_classif.shape)

    blastx_classif = merged_bx[~merged_bx['IndividualID'].isin(blastn_classif['IndividualID'])]
    print(blastx_classif.shape)

    merged = pd.concat([blastn_classif, blastx_classif])
    print(merged.shape)

    merged_missed = merged[merged['lineage'].isnull()]
    merged_missed.to_csv(out_merged_missed, sep='\t', index=False)

    # drop missed from merged with copy
    merged_classif = merged.copy()[merged['lineage'].notnull()]
#    merged_classif = merged[merged['lineage'].notnull()]
#    merged_classif.to_csv(out_besthit_classif, sep='\t', index=False)

    lineage_columns = ['superkingdom', 'realm', 'subrealm', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'suborder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'name']
    merged_classif[lineage_columns] = merged_classif['lineage'].str.split(',', expand=True) # split lineage column into separate columns
    merged_classif['name'] = merged_classif['name'].str.replace('_', ',') #put back any commas in the name
    merged_classif.to_csv(out_besthit_classif, sep='\t', index=False)
    
    # produce a besthit with only a select number of columns ('clean')
    vmr_cols_to_keep = ['order', 'suborder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'name']
    renamed_vmr = ['blastvmr_' + x for x in vmr_cols_to_keep]
    renameD = dict(zip(vmr_cols_to_keep, renamed_vmr))
    cols_to_keep = existing_spillover_columns + vmr_cols_to_keep
    shortmerged = merged_classif[cols_to_keep]
    shortmerged = shortmerged.rename(columns=renameD)
    shortmerged.to_csv(out_besthit_clean, sep='\t', index=False)


def cmdline(sys_args):
    p = argparse.ArgumentParser()
    p.add_argument('--spillover-csv', help='spillover info file', default='inputs/2023-03-27_spillover_accession-numers.csv')
    p.add_argument('-n', '--blastn', help='Blastn file', default='output.spillover-blast/blastn/spillover-x-vmr_MSL38_v1.blastn.best.tsv')
    p.add_argument('-x', '--blastx', help='Blastx file', default='output.spillover-blast/diamond-blastx/spillover-x-vmr_MSL38_v1.diamond-blastx.best.tsv')
    p.add_argument('-o', '--output-base' ,help='Output base', required=True)
    args = p.parse_args()
    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
