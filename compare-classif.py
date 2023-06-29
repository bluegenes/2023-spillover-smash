import pandas as pd
gather_file = 'spillover.dna-k31.gather-classifications.csv'
blast_file = 'sfv-x-ictv-spumavirus.best.tsv'
gDF = pd.read_csv(gather_file, sep = ',')
bDF = pd.read_csv(blast_file, sep = '\t')
bDF.rename(columns = {'qseqid': 'ident', 'taxonomy': 'blast_lineage'}, inplace=True)
bDF['blast_lineage'] = bDF['blast_lineage'].str.rsplit(';', n=1, expand=True)[0]
blastn = bDF.loc[bDF['source'] == 'blastn']
# split virus name from blast (we don't have in gather)


gDF['lineage'] = 'Viruses;' + gDF['lineage']
gDF.rename(columns = {'lineage': 'gather-k31-lineage'}, inplace=True)
gDF['ident'] = gDF['query_name'].str.split(' ', expand=True)[0]

cDF = blastn.merge(gDF, on='ident', how='left')

cDF.to_csv('spillover.classif-combined.csv', index=False)

# count number with same lineage
cDF[['blast_lineage', 'gather-k31-lineage']]
count = sum(cDF['blast_lineage'] == cDF['gather-k31-lineage'])
print(count)
cDF['blast_lineage'][0]
