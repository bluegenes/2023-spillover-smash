import sourmash
import pandas as pd
import numpy as np 

from sourmash.tax.tax_utils import ICTVRankLineageInfo

# read in spillover csv
sp_csv = '../inputs/2023-03-27_spillover_accession-numers.csv'
sDF = pd.read_csv(sp_csv)
#sDF

compare_csv =  'spillover.classif-combined.csv'
cDF = pd.read_csv(compare_csv)

# make ICTVLineageInfo object for each lineage
def make_ICTVLineageInfo(row, column='blast_lineage'):
    lin = row[column]
    if not isinstance(lin, str): # if nan
        return ICTVRankLineageInfo()
    lin_str = lin.split(';', 1)[1] # drop 'Viruses;'
    lineage =  ICTVRankLineageInfo(lineage_str = lin_str)
    return lineage


# make sourmash lineages
cDF['blast_linfo'] = cDF.apply(make_ICTVLineageInfo, axis=1, column='blast_lineage')
cDF['gather-k31-linfo'] = cDF.apply(make_ICTVLineageInfo, axis=1, column='gather-k31-lineage')

# find LCA
def find_lca(row, col1='blast_linfo',col2='gather-k31-linfo'):
    lin1 = row[col1]
    lin2 = row[col2]
    lca_info = lin1.find_lca(lin2)
    if lca_info is None:
        row['lca'] = ''
        row['lca_rank'] = ''
    else:
        row['lca'] = lca_info.display_lineage()
        row['lca_rank'] = lca_info.lowest_rank
    return row

cDF = cDF.apply(find_lca, axis=1, col1='blast_linfo', col2='gather-k31-linfo')

# count number with each lca_rank
count = cDF['lca_rank'].value_counts()
print(count)

# get average percent identity and mismatch by lca_rank
cDF.groupby(['lca_rank']).agg({'pident': ['mean', 'std']})
# modify the above to average both pident and mismatch
cDF.groupby(['lca_rank']).agg({'pident': ['mean', 'median', 'std'], 
                               'mismatch': ['mean', 'median', 'std'],
                            'f_weighted_at_rank': ['mean', 'median', 'std']})


cDF.groupby(['lca_rank']).agg({'pident': ['mean', 'std']})
# mismatch


write_cDF = cDF.drop(columns=['blast_linfo', 'gather-k31-linfo', 'source'], inplace=True)
write_cDF.to_csv('spillover.classif-combined.assess.csv', index=False)

# extract lineage to columns
def lin_to_cols(row, column='blast_linfo', ranks_to_keep = ['name','species', 'genus', 'family']):
    prefix = 'blast_vmr_'
    if 'gather' in column:
        prefix = 'gather_vmr_'
    lin = row[column]
    lin_dict = {prefix + a.rank: a.name for a in lin.lineage}
    for k, v in lin_dict.items():
        if k.split(prefix)[1] in ranks_to_keep:
            row[k] = v
    return row

blast_columns = ['ident', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
       'qend', 'sstart', 'send', 'evalue', 'bitscore', 'blast_lineage', 'blast_linfo']
blastDF = cDF.copy()[blast_columns]
blastDF.rename(columns= {'ident': 'qseqid'}, inplace=True) # restore qseqid name
blastDF = blastDF.apply(lin_to_cols, axis=1, column='blast_linfo')

blast_merged = sDF.merge(blastDF, right_on='qseqid', left_on='AccessionNumber', how='right')

blast_merged_write = blast_merged.drop(columns=['blast_linfo'])
blast_merged_write.to_csv('spillover.classif-combined.assess.blast-only.csv', index=False)


# Baboon comparisons
baboon_csv = 'baboon-colony-africa-taiwan-chk.csv' 
bDF = pd.read_csv(baboon_csv)

keep_columns = ['IndividualID', 'AccessionNumber', 'Virus', 'VirusGenus',
       'VirusSpecies', 'VirusFamily', 'HostSpecies', 'HostGenus', 'HostFamily',
        'Virus_vmr', 'Virus_orig', 'bp', 'HostSpecies_common_name']

bDF = bDF[keep_columns]

# merge with blast, gather results
merged = bDF.merge(cDF, left_on='AccessionNumber', right_on='ident', how='left')

merged_plus = merged.apply(lin_to_cols, axis=1, column='blast_linfo')
merged_plus = merged_plus.apply(lin_to_cols, axis=1, column='gather-k31-linfo')

merged_plus[['AccessionNumber','blast_vmr_species', 'gather_vmr_species', 'HostSpecies_common_name']]

merged_plus_write = merged_plus.drop(columns=['blast_linfo', 'gather-k31-linfo', 'source'])
merged_plus_write.to_csv('baboon-colony-africa-taiwan-chk.assess.csv', index=False)

blast_only_columns_reordered =  ['IndividualID', 'AccessionNumber', 'Virus', 'VirusGenus',
       'VirusSpecies', 'VirusFamily', 'HostSpecies', 'HostGenus', 'HostFamily',
       'Virus_orig', 'bp', 'HostSpecies_common_name',
       'blast_vmr_name', 'blast_vmr_species',
       'blast_vmr_genus', 'blast_vmr_family',
       'sseqid', 'pident',  'length', 'mismatch', 
       'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
       'bitscore', 'blast_lineage']

merged_blast_only = merged_plus_write[blast_only_columns_reordered]
merged_blast_only.to_csv('baboon-colony-africa-taiwan-chk.assess.blast-only.csv', index=False)


# # extract blast info to separate columns
# def lineages_to_select_columns(row, column='blast_linfo'):
#     lin = row[column]
#     new_name_col = 'blast_vmr_name'
#     new_species_col = 'blast_vmr_species'
#     new_genus_col = 'blast_vmr_genus'
#     new_family_col = 'blast_vmr_family'
#     if 'gather' in column:
#         new_name_col = 'gather_vmr_virus_name'
#         new_species_col = 'gather_vmr_virus_species'
#         new_genus_col = 'gather_vmr_virus_genus'
#         new_family_col = 'gather_vmr_virus_family'
#     if not isinstance(lin, ICTVRankLineageInfo): # if nan
#         row[new_name_col] = ''
#         row[new_species_col] = ''
#         row[new_genus_col] = ''
#         return row
#     row[new_name_col] = lin.name_at_rank('name')
#     row[new_species_col] = lin.name_at_rank('species')
#     row[new_genus_col] = lin.name_at_rank('genus')
#     row[new_family_col] = lin.name_at_rank('family')
#     return row