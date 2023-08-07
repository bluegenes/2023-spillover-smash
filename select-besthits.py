import os, sys
import argparse
import re
import pandas as pd

def add_taxonomy(row, gene2lineage):
    match_acc = row['sseqid']
    # blastn results have a prefix in front of accession, e.g. 'ref|NC_001422.1|'
    if any(x in match_acc for x in ['ref|', 'gb|', 'dbj|', 'emb|']):
#    if '|' in match_acc: # some blastn results have 'ref|' in front of accession
        pattern = r'\|([^|]+)\|'
        match1 = re.search(pattern, match_acc)
        if match1:
            match_acc = match1.group(1)
        else:
            raise ValueError('Could not parse ref| accession from blastn result.')
    match_acc = match_acc.rsplit('.')[0]
    if match_acc not in gene2lineage.keys():
        if args.force:
            print(f"WARNING: {match_acc} lineage not found.")
            row['lineage'] = 'NA'
        else:
            raise ValueError(f"ERROR: {match_acc} lineage not found.")
    else:
        row['lineage'] = gene2lineage[match_acc]
    return row

def expand_gene_level_accs(row):
    gene_accs, prot_accs = [], []
    if pd.notnull(row['gene_accs']):
        gene_accs = row['gene_accs'].split(',') if ',' in row['gene_accs'] else [row['gene_accs']]
        gene_accs = [str(x).split('.')[0] for x in gene_accs]
    if pd.notnull(row['protein_accs']):
        protein_accs = row['protein_accs'].split(',') if ',' in row['protein_accs'] else [row['protein_accs']]
        prot_accs = [str(x).split('.')[0] for x in protein_accs]
    row['gene_accessions'] = gene_accs
    row['protein_accessions'] = prot_accs
    return row

def select_best_hits(df):
    # select top bitscore as best hit
    return df.dropna(subset=['bitscore']).groupby(['qseqid', 'source']).apply(lambda x: x.loc[x['bitscore'].idxmax()]).reset_index(drop=True)

def main(args):
    # make output filenames
    out_all = args.output_base + '.all.tsv'
    out_best = args.output_base + '.best.tsv'

    blast_file = args.blast
    blast_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(' ')
    blast_df = pd.read_csv(blast_file, sep='\t', names=blast_columns)
    blast_df['source'] = args.searchtype

    taxDF = pd.read_csv(args.tax_file, sep='\t', dtype={'gene_accs': str, 'protein_accs': str, 'lineage': str})
    taxDF = taxDF.apply(expand_gene_level_accs, axis=1)
    gene_exploded = taxDF.explode(['gene_accessions'])
    prot_exploded = taxDF.explode(['protein_accessions'])
    gene2lineage = dict(zip(gene_exploded['gene_accessions'], gene_exploded['lineage']))
    prot2lineage = dict(zip(prot_exploded['protein_accessions'], prot_exploded['lineage']))
    gene2lineage.update(prot2lineage)

    blast_df = blast_df.apply(add_taxonomy, args=(gene2lineage,), axis=1)

    # merge with spillover info
    sDF = pd.read_csv(args.spillover_csv, index_col=0)
    blast_merged = sDF.merge(blast_df, right_on='qseqid', left_on='AccessionNumber', how='left')
    blast_merged.to_csv(out_all, sep='\t', index=False)

    # select best hits
    best = select_best_hits(blast_merged)
    best.to_csv(out_best, sep='\t', index=False)

def cmdline(sys_args):
    p = argparse.ArgumentParser()
    p.add_argument('--blast', help='tab-separated blast output (outfmt 6) or diamond blast output', required=True)
    p.add_argument('-t', '--tax-file', help='taxonomy TSV file', default="output.vmr/vmr_MSL38_v1.taxonomy.tsv") 
    p.add_argument('--spillover-csv', help='spillover info file', default='inputs/2023-03-27_spillover_accession-numers.csv')
    p.add_argument('--searchtype', help='search type')
    p.add_argument('--force', help='force continue if lineage not found', action='store_true')
    p.add_argument('-o', '--output-base' ,help='Output base', default=None)
    args = p.parse_args()
    
    if not args.searchtype:
        # look for 'blastn' or 'diamond-blastx' in blast filename
        if 'blastn' in args.blast:
            args.searchtype = 'blastn'
        elif 'diamond-blastx' in args.blast:
            args.searchtype = 'diamond-blastx'
        elif 'blastx' in args.blast:
            args.searchtype = 'blastx'
    if not args.output_base:
        args.output_base = args.blast.rsplit(f'.tsv')[0]
    return main(args)

# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
