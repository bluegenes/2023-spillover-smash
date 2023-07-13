import pandas as pd
import csv
import numpy as np

db_dir = "output.vmr"
out_dir = "output.spillover-blast"
logdir = os.path.join(out_dir, 'logs')

# get genome filenames
fromfile = 'output.spillover/spillover.fromfile.csv'
ff = pd.read_csv(fromfile)
ff['ident'] = ff['name'].str.split(' ', expand=True)[0]
# drop rows if ident in null list
null_list = ["", np.nan, 'Unknown']
ff = ff[~ff['ident'].isin(null_list)]
print(ff.shape)

dna_acc = ff["ident"].tolist()
prot_ff = ff[~ff['protein_filename'].isna()]#['ident']
ff.set_index('ident', inplace=True)

# full set
basename = 'spillover'
db_basename = 'vmr_MSL38_v1'
ACCS = dna_acc
PROT_ACCS = prot_ff['ident'].tolist()

# test simian foamy virus
#basename="sfv"
#db_basename = "ictv-spumavirus"
#sm_file = 'inputs/simian-foamy.spillover.csv'
#sm = pd.read_csv(sm_file)
#ACCS = sm['AccessionNumber']

# test set lassa
#basename="lassa"
#db_basename= 'vmr_MSL38_v1.lassa'
#sm_file  = 'inputs/lassa.spillover.csv'
#sm = pd.read_csv(sm_file)
#ACCS = sm['AccessionNumber']

# get fastas from ff file for combine_fasta
dna_fastas = ff["genome_filename"].tolist()
prot_fastas = prot_ff["protein_filename"].tolist()
fastaD = {'dna': dna_fastas, 'protein': prot_fastas}

# # dna taxonomy
# tax_file = 'inputs/VMR_MSL38_v1.acc.csv'
# ictv_tax = pd.read_csv(tax_file)
# #df['taxonomy'] = df[columns_to_select].apply(lambda x: ';'.join(x.dropna().astype(str)), axis=1)
# ictv_tax['taxonomy'] = 'Viruses;' + ictv_tax[columns_to_select].fillna('').apply(lambda x: ';'.join(x.astype(str)), axis=1)

# ictvD = dict(zip(ictv_tax['Virus REFSEQ accession'], ictv_tax['taxonomy']))
# ictvD_gb = dict(zip(ictv_tax['Virus GENBANK accession'], ictv_tax['taxonomy']))
# ictvD.update(ictvD_gb) # add gb mappings

# use database taxonomy files instead
tax_file = os.path.join(db_dir, f"{db_basename}.taxonomy.csv")
#prot_tax_file = os.path.join(db_dir, f"{db_basename}.protein-taxonomy.csv")

# to do: use protein?
#PROT_ACCS = [x for x in ACCS if x in PROT_EXISTS]
#n_missing = len(ACCS) - len(PROT_ACCS)
#print(f"missing {n_missing} protein accessions")

onstart:
    print("------------------------------")
    print("blast workflow")
    print("------------------------------")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")

rule all:
    input:
#        os.path.join(out_dir, "blastn-combined", f"{basename}-x-{db_basename}.blastn.tsv"),
#        os.path.join(out_dir, 'diamond-combined', f"{basename}-x-{db_basename}.diamond-blastx.tsv"),
        expand(os.path.join(out_dir, f"{basename}-x-{db_basename}.{{searchtype}}.{{end}}.tsv"), searchtype = ['blastn', 'diamond-blastx'], end = ['best', 'all']), # best.info
        os.path.join(out_dir, f"{basename}-x-{db_basename}.blast-classif.tsv"),
        os.path.join(out_dir, f"{basename}-x-{db_basename}.missed.tsv"),


rule combine_fasta:
    input: 
         fastas= lambda w: fastaD[w.moltype],
    output: os.path.join(out_dir,  "combined", "{basename}.{moltype}.fa.gz")
    params:
    log:  os.path.join(logdir, "combined", "{basename}.{moltype}.log")
    benchmark:  os.path.join(logdir, "combined", "{basename}.{moltype}.benchmark")
    shell:
        """
        zcat {input.fastas} | gzip > {output} 2> {log}
        """


rule blastn:
    input:
        query =  os.path.join(out_dir, "combined","{basename}.dna.fa.gz"),
        index = os.path.join(db_dir, "blast", "{db_basename}.dna.index.nhr"),
    output:
        results = os.path.join(out_dir, "blastn", "{basename}-x-{db_basename}.blastn.tsv")
    params:
        index_basename = os.path.join(db_dir, "blast", "{db_basename}.dna.index")
    log: os.path.join(logdir, 'blastn', "{basename}-x-{db_basename}.blastn.log")
    benchmark: os.path.join(logdir, 'blastn', "{basename}-x-{db_basename}.blastn.benchmark")
    conda: "conf/env/blast.yml"
    threads: 30
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        time=10000,
        partition="bmh",
    shell:
        """
        zcat {input.query} | blastn -query - -db {params.index_basename} \
               -out {output.results} -num_threads {threads} \
               -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
               2> {log}
        """

rule diamond_blastx:
    input:
        query = os.path.join(out_dir, "combined", "{basename}.dna.fa.gz"), 
        index = os.path.join(db_dir, "diamond", "{db_basename}.protein.fa.gz" + ".dmnd"),
    output:
        results = os.path.join(out_dir, 'diamond-blastx', "{basename}-x-{db_basename}.diamond-blastx.tsv"),
    conda: "conf/env/diamond.yml"
    log: os.path.join(logdir, 'diamond-blastx', "{basename}-x-{db_basename}.diamond-blastx.log")
    benchmark: os.path.join(logdir, 'diamond-blastx', "{basename}-x-{db_basename}.diamond-blastx.benchmark")
    threads: 30
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        time=6000,
        partition="high2",
    shell:
        """
        diamond blastx --query {input.query} --db {input.index} \
                       --out {output.results} --sensitive \
                       --threads {threads} 2> {log}
        """

localrules: select_besthits
# select besthits and add taxonomy
rule select_besthits:
    input:
        blast = os.path.join(out_dir, "{searchtype}", "{basename}-x-{db_basename}.{searchtype}.tsv"),
        spillover_info = 'inputs/2023-03-27_spillover_accession-numers.csv',
    output:
        all = os.path.join(out_dir, "{basename}-x-{db_basename}.{searchtype}.all.tsv"),
        best = os.path.join(out_dir, "{basename}-x-{db_basename}.{searchtype}.best.tsv"),
#        all_info = os.path.join(out_dir, "{basename}-x-{db_basename}.{searchtype}.all.info.tsv"),
 #       best_info = os.path.join(out_dir, "{basename}-x-{db_basename}.{searchtype}.best.tsv"),
    run:
        import re
        def select_best_hits(df):
            return df.dropna(subset=['bitscore']).groupby(['qseqid', 'source']).apply(lambda x: x.loc[x['bitscore'].idxmax()]).reset_index(drop=True)
        
        def add_taxonomy(row):
            match_acc = row['sseqid']
            if '|' in match_acc:
                pattern = r'\|([^|]+)\|' # match the first group of characters between two pipes
                match1 = re.search(pattern, match_acc).group(1)
                match_acc = match1
            match_acc = match_acc.rsplit('.')[0]
            if match_acc not in gene2lineage.keys():
                print(f"WARNING: {match_acc} lineage not found.")
                row['lineage'] = 'NA'
            else:
                row['lineage'] = gene2lineage[match_acc] 
            return row
        
        def expand_gene_level_accs(row):
           # convert ';'-separated list of gene accessions to list
            gene_accs, prot_accs = [], []
            if pd.notnull(row['gene_accs']):
                gene_accs = row['gene_accs'].split(';') if ';' in row['gene_accs'] else [row['gene_accs']]
                gene_accs = [str(x).split('.')[0] for x in gene_accs] # remove '.1' from gene accessions
            else:
                # these should reflect 'suppressed records' where we don't have a genome file
                print(f"WARNING: gene {row['gene_accs']} is not a string")
            if pd.notnull(row['protein_accs']):
                protein_accs = row['protein_accs'].split(';') if ';' in row['protein_accs'] else [row['protein_accs']]
                prot_accs = [str(x).split('.')[0] for x in protein_accs] # remove '.1' from protein accessions
            else:
                # these should represent anywhere we don't have a proteome file (suppressed records and others)
                print(f"WARNING: protein {row['protein_accs']} is not a string")
            row['gene_accessions'] = gene_accs
            row['protein_accessions'] = prot_accs
            return row

        blast_file = input.blast
        # Read the blast file into a DataFrame
        blast_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(' ')
        blast_df = pd.read_csv(blast_file, sep='\t', names=blast_columns)
        
        # Add an additional column to indicate the source (blastn or blastx)
        blast_df['source'] = wildcards.searchtype

        # create dictionary of gene-level accessions to genome accession
        taxDF = pd.read_csv(tax_file, dtype={'gene_accs': str, 'protein_accs': str, 'lineage': str})
        taxDF = taxDF.apply(expand_gene_level_accs, axis=1)
        gene_exploded = taxDF.explode(['gene_accessions'])
        prot_exploded = taxDF.explode(['protein_accessions'])
        gene2lineage = dict(zip(gene_exploded['gene_accessions'], gene_exploded['lineage']))
        prot2lineage= dict(zip(prot_exploded['protein_accessions'], prot_exploded['lineage']))
        gene2lineage.update(prot2lineage) # add protein mappings

        # use apply with add_taxonomy to add lineage info to blast_df
        blast_df = blast_df.apply(add_taxonomy, axis=1)
        # blast_df.to_csv(output.all, sep='\t', index=False)


        # add spillover info
        sDF = pd.read_csv(str(input.spillover_info), index_col= 0)
        blast_merged = sDF.merge(blast_df, right_on='qseqid', left_on='AccessionNumber', how='left') # 'left' to keep all spillover accessions

        blast_merged.to_csv(output.all, sep='\t', index=False)
        
        # select best per df
        best = select_best_hits(blast_df)
        best.to_csv(output.best, sep='\t', index=False)
        # add spillover info
        #select columns
#       ranks_to_keep = ['name', 'species', 'genus', 'family']
#       merged_write = merged.drop(columns=['blast_linfo'])
#       merged_write.to_csv(args.output, sep = '\t', index=False)



rule combine_best:
    input:
        blastn = os.path.join(out_dir, "{basename}-x-{db_basename}.blastn.best.tsv"),
        blastx = os.path.join(out_dir, "{basename}-x-{db_basename}.diamond-blastx.best.tsv"),
    output:
        besthit_classif = os.path.join(out_dir, "{basename}-x-{db_basename}.blast-classif.tsv"),
        missed = os.path.join(out_dir, "{basename}-x-{db_basename}.missed.tsv"),
    # conda: 'conf/env/reports.yml'
    run:
        # read blastn as dataframe
        blastn_df = pd.read_csv(input.blastn, sep='\t')
        # read blastx as dataframe
        blastx_df = pd.read_csv(input.blastx, sep='\t')
        # if lineage is null in blastn, replace with lineage from blastx
        blastn_df.loc[blastn_df['lineage'].isnull(), 'lineage'] = blastx_df['lineage']
        # if source is null but lineage has been filled, set source to blastx
        blastn_df.loc[blastn_df['source'].isnull() & blastn_df['lineage'].notnull(), 'source'] = 'blastx'
        # if lineage is still null, write to missed csv
        blastn_missed = blastn_df[blastn_df['lineage'].isnull()]
        blastn_missed.to_csv(output.missed, sep='\t', index=False)

        # write to blast-classif
        blastn_df.to_csv(output.besthit_classif, sep='\t', index=False)
