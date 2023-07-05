import pandas as pd

vmr_file = 'inputs/VMR_MSL38_v1.acc.csv'
tax_file = 'inputs/VMR_MSL38_v1.taxonomy.csv'
# tax_file = 'VMR_21-221122_MSL37.taxonomy.csv'
vmr = pd.read_csv(vmr_file)
vmr = vmr.rename(columns={'GenBank Assembly ID':'ident', 'Virus name(s)': 'name', 'Exemplar or additional isolate': 'exemplar_or_additional'})
print(vmr.shape)
vmr = vmr.dropna(subset=['ident'])
vmr = vmr.drop_duplicates(subset=['ident']) #only keep first
print(vmr.shape)
vmr.set_index('ident', inplace=True)
vmr['superkingdom'] = 'Viruses'
tax_columns = ['superkingdom', 'Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', \
               'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', \
               'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species', \
               'exemplar_or_additional', 'name']

tax_info = vmr[tax_columns]
tax_info = tax_info.rename(columns=str.lower) # lowercase all names
tax_info.to_csv(tax_file)
