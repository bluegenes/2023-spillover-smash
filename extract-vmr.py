# Script to get long form VMR csvs with segment-specific accessions

import pandas as pd
import numpy as np

# VMR Reference file, v37
vmr_file = "inputs/VMR_21-221122_MSL37.xlsx"
vmr = pd.read_excel(vmr_file, sheet_name="VMRb37")

# for segmented viruses, there's more than one accession. How do we want to handle this from taxonomy point of view?
# signature names: "{accession} {name} {segment}"

# other important info:
# column: "Exemplar or additional isolate"? 'E' = 'exemplar', 'A' = 'additional'
# Taxonomic columns:
# Realm, Subrealm, Kingdom, Subkingdom, Phylum, Subphylum, Class, Subclass, Order, Suborder, Family, Subfamily, Genus, Subgenus, Species, Virus name(s)

# explode genbank/refseq info columns. keep these separate for now. Think about making genbank: refseq map file?
gb = vmr.assign(genbank=vmr['Virus GENBANK accession'].str.split(';')).explode(['genbank'])
rf = vmr.assign(refseq=vmr['Virus REFSEQ accession'].str.split(';')).explode(['refseq'])
gb.reset_index(drop=True, inplace=True)
rf.reset_index(drop=True, inplace=True)

#some entries have segnments, some do not. segments have format name: acc
# use str.extract to extract the values before and after the ':':
gb[['segment', 'accession']] = gb['genbank'].str.extract('(.*:)?\s?(.*)')
rf[['segment', 'accession']] = rf['refseq'].str.extract('(.*:)?\s?(.*)')
# now strip colons and whitespace and replace nan with ''
gb['segment'] = gb['segment'].str.rstrip(':').str.strip().fillna("")
rf['segment'] = rf['segment'].str.rstrip(':').str.strip().fillna("")
# these aren't really needed, but just in case, strip whitespace from acc col too
gb['accession'] = gb['accession'].str.strip().fillna("")
rf['accession'] = rf['accession'].str.strip().fillna("")
# ok build signature names from this info
gb['signame'] = gb['accession'] + ' ' + gb['Virus name(s)'] + ' ' + gb['segment']
rf['signame'] = rf['accession'] + ' ' + rf['Virus name(s)'] + ' ' + rf['segment']

#strip any trailing spaces introduced by building names without segments
gb['signame'] = gb['signame'].str.strip()
rf['signame'] = rf['signame'].str.strip()

# write new full dataframes
gb.to_csv("inputs/VMR_21-221122_MSL37.genbank.csv", index=False)
rf.to_csv("inputs/VMR_21-221122_MSL37.refseq.csv", index=False)

# build taxonomy file
tax_columns = ['accession', 'Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species', 'Virus name(s)', 'Virus name abbreviation(s)', 'Virus isolate designation', 'signame', 'Exemplar or additional isolate']
tax_gb = gb.copy()[tax_columns]
tax_rf = rf.copy()[tax_columns]

# rename exemplar column
renameD = {'Exemplar or additional isolate': 'exemplar',
            'Virus name(s)': 'name',
            'Virus name abbreviation(s)': 'abbreviation',
            'Virus isolate designation': 'isolate_designation' }
tax_gb.rename(columns=renameD, inplace=True)
tax_rf.rename(columns=renameD, inplace=True)

# use T/F instead of A/E in exemplar column
tax_gb['exemplar'] = tax_gb['exemplar'].replace({'E': True, 'A': False})
tax_rf['exemplar'] = tax_gb['exemplar'].replace({'E': True, 'A': False})

# some entries don't have accessions! We can't use these -- let's drop from taxonomy info
# first turn the "" back into np.nan
tax_gb['accession'].replace('', np.nan, inplace=True)
tax_rf['accession'].replace('', np.nan, inplace=True)
# now drop these rows
tax_gb = tax_gb.dropna(subset=['accession'])
tax_rf = tax_rf.dropna(subset=['accession'])

# write new dataframes
tax_gb.to_csv("inputs/VMR_21-221122_MSL37.genbank.taxonomy.csv", index=False)
tax_rf.to_csv("inputs/VMR_21-221122_MSL37.refseq.taxonomy.csv", index=False)
