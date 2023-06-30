import pandas as pd
import csv

db_dir = "output.vmr"
out_dir = "output.spillover-blast"
logdir = os.path.join(out_dir, 'logs')



# get genome filenames
fromfile = 'output.spillover/spillover.fromfile.csv'
ff = pd.read_csv(fromfile)
ff['ident'] = ff['name'].str.split(' ', expand=True)[0]
dna_acc = sp["ident"][sp["genome_filename"].notnull()].tolist()
PROT_EXISTS = ff[~ff['protein_filename'].isna()]['ident']
ff.set_index('ident', inplace=True)

# full set
basename = 'spillover'
db_basename = 'ictv'
ACCS = dna_acc

# test simian foamy virus
basename="sfv"
db_basename = "ictv-spumavirus"
sm_file = 'inputs/simian-foamy.spillover.csv'
sm = pd.read_csv(sm_file)
ACCS = sm['AccessionNumber']


# ictv taxonomy
tax_file = 'inputs/VMR_MSL38_v1.acc.csv'
ictv_tax = pd.read_csv(tax_file)
columns_to_select = ['Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Subclass',
                     'Order', 'Suborder', 'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species', 'Virus name(s)']
#df['taxonomy'] = df[columns_to_select].apply(lambda x: ';'.join(x.dropna().astype(str)), axis=1)
ictv_tax['taxonomy'] = 'Viruses;' + ictv_tax[columns_to_select].fillna('').apply(lambda x: ';'.join(x.astype(str)), axis=1)

ictvD = dict(zip(ictv_tax['Virus REFSEQ accession'], ictv_tax['taxonomy']))
ictvD_gb = dict(zip(ictv_tax['Virus GENBANK accession'], ictv_tax['taxonomy']))
ictvD.update(ictvD_gb) # add gb mappings



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
        expand(os.path.join(out_dir, 'blast', "{acc}-x-{db_basename}.blast.txt"), acc = ACCS ,db_basename=db_basename),
        expand(os.path.join(out_dir, 'diamond', "{acc}-x-{db_basename}.diamond-blastx.txt"), acc = ACCS, db_basename=db_basename), 
        os.path.join(out_dir, f"{basename}-x-{db_basename}.blast.tsv"),
        os.path.join(out_dir, f"{basename}-x-{db_basename}.diamond-blastx.tsv"),
        os.path.join(out_dir, f"{basename}-x-{db_basename}.best.tsv"),
#        os.path.join(out_dir, f"{basename}-x-{db_basename}.blast.best.tsv"),

rule blastn:
    input:
        query = lambda w: ff.at[w.acc, 'genome_filename'],
        index = os.path.join(db_dir, "blast", "{db_basename}.dna.index.nhr")
    output:
        results = os.path.join(out_dir, 'blast', "{acc}-x-{db_basename}.blast.txt")
    params:
        index_basename = os.path.join(db_dir, "blast", "{db_basename}.dna.index")
    log: os.path.join(logdir, 'blastn', "{acc}-x-{db_basename}.blastn.log")
    benchmark: os.path.join(logdir, 'blastn', "{acc}-x-{db_basename}.blastn.benchmark")
    conda: "conf/env/blast.yml"
    threads: 4
    shell:
        """
        zcat {input.query} | blastn -query - -db {params.index_basename} \
               -out {output.results} -num_threads {threads} \
               -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
               2> {log}
        """


rule diamond_blastx:
    input:
        query = lambda w: ff.at[w.acc, 'genome_filename'],
        index = os.path.join(db_dir, "diamond", "{db_basename}.protein.fa.gz" + ".dmnd"),
    output:
        results = os.path.join(out_dir, 'diamond', "{acc}-x-{db_basename}.diamond-blastx.txt"),
    conda: "conf/env/diamond.yml"
    log: os.path.join(logdir, 'diamond', "{acc}-x-{db_basename}.diamond-blastx.log")
    benchmark: os.path.join(logdir, 'diamond', "{acc}-x-{db_basename}.diamond-blastx.benchmark")
    threads: 4
    shell:
        """
        diamond blastx --query {input.query} --db {input.index} \
                       --out {output.results} --sensitive \
                       --threads {threads} 2> {log}
        """


rule diamond_blastp:
    input:
        query = lambda w: ff.at[w.acc, 'protein_filename'],
        index = os.path.join(db_dir, "diamond", "{db_basename}.protein.fa.gz" + ".dmnd"),
    output:
        results = os.path.join(out_dir, 'diamond', "{acc}-x-{db_basename}.diamond-blastp.txt"),
    conda: "conf/env/diamond.yml"
    log: os.path.join(logdir, 'diamond', "{acc}-x-{db_basename}.diamond-blastp.log")
    benchmark: os.path.join(logdir, 'diamond', "{acc}-x-{db_basename}.diamond-blastp.benchmark")
    threads: 4
    shell:
        """
        diamond blastp --query {input.query} --db {input.index} \
                       --out {output.results} --quiet \
                       --sensitive --threads {threads} 2> {log}
        """

rule aggregate_blastn:
    input:
        cl = lambda w: expand(os.path.join(out_dir, 'blast', "{acc}-x-{db_basename}.blast.txt"), acc = ACCS, db_basename=db_basename),
    output:
        os.path.join(out_dir, "{basename}-x-{db_basename}.blast.tsv")
    run:
        with open(str(output), 'w') as outF:
            w = csv.writer(outF, delimiter='\t')
            header = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore taxonomy".split(' ')
            w.writerow(header)
            empties = []
            for inF in input:
                inF = str(inF)
                if os.stat(inF).st_size == 0: # file is empty
                    acc = os.path.basename(inF).split('-', 1)[0]
                    empties.append([acc])
                with open(inF, "r", newline="") as cl:
                    r = csv.reader(cl,delimiter='\t')
                    for row in r:
                        if '|' in row[1]:
                            import re
                            pattern = r'\|([^|]+)\|'
                            match1 = re.search(pattern, row[1]).group(1)
                            row[1] = match1
                        match_acc = row[1]
                        match_acc = match_acc.rsplit('.')[0]
                        taxonomy = ictvD.get(match_acc, '')
                        #taxonomy = ictv_tax.at[match_acc, "taxonomy"]
                        row.append(taxonomy)
                        w.writerow(row)
            for e in empties:
                w.writerow(e)


rule aggregate_blastx:
    input:
        cl = lambda w: expand(os.path.join(out_dir, 'diamond', "{acc}-x-{db_basename}.diamond-blastx.txt"), acc = ACCS, db_basename=db_basename),
    output:
        os.path.join(out_dir, "{basename}-x-{db_basename}.diamond-blastx.tsv"),
    run:
        with open(str(output), 'w') as outF:
            w = csv.writer(outF, delimiter='\t')
            header = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore taxonomy".split(' ')
            w.writerow(header)
            empties = []
            for inF in input:
                inF = str(inF)
                if os.stat(inF).st_size == 0: # file is empty
                    acc = os.path.basename(inF).split('-', 1)[0]
                    empties.append([acc])
                with open(inF, "r", newline="") as cl:
                    r = csv.reader(cl,delimiter='\t')
                    for row in r:
                        match_acc = row[1]
                        match_acc = match_acc.rsplit('.')[0]
                        #taxonomy = ictv_tax.at[match_acc, "taxonomy"]
                        taxonomy = ictvD.get(match_acc, '')
                     #   taxonomy = ictv_tax.at[row[1], "taxonomy"]
                     #   row.append(taxonomy)
                        w.writerow(row)
            for e in empties:
                w.writerow(e)

rule select_besthits:
    input:
        blastx = os.path.join(out_dir, "{basename}-x-{db_basename}.diamond-blastx.tsv"),
        blastn = os.path.join(out_dir, "{basename}-x-{db_basename}.blast.tsv"),
    output:
        os.path.join(out_dir, "{basename}-x-{db_basename}.best.tsv"),
    run:
        def select_best_hits(df):
            return df.dropna(subset=['bitscore']).groupby(['qseqid', 'source']).apply(lambda x: x.loc[x['bitscore'].idxmax()]).reset_index(drop=True)
        
        blastn_file = input.blastn
        blastx_file = input.blastx
        output_file = output[0]
        
        # Read the blastn and blastx files into DataFrames
        blastn_df = pd.read_csv(blastn_file, sep='\t')
        blastx_df = pd.read_csv(blastx_file, sep='\t')
        
        # Add an additional column to indicate the source (blastn or blastx)
        blastn_df['source'] = 'blastn'
        blastx_df['source'] = 'blastx'
        
        # Combine the blastn and blastx DataFrames
        combined_df = pd.concat([blastn_df, blastx_df])
        
        # Select rows with the maximum bit score per qseqid
        best_hits_df = select_best_hits(combined_df)
        #best_hits_df = select_best_hits(blastn_df)
        
        # Write the best hits to the output file
        best_hits_df.to_csv(output_file, sep='\t', index=False)
