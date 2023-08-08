import pandas as pd
import csv
import numpy as np

# to do: set these from config file
db_dir = "output.vmr"
out_dir = "output.spillover-blast"
logdir = os.path.join(out_dir, 'logs')

#spillover_csv = 'inputs/2023-03-27_spillover_accession-numers.csv'
spillover_csv = 'inputs/2023-08-07_spillover_accession-numbers_info.csv'

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
basename = '2023-08-07-spillover'
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

# use database taxonomy files instead
tax_file = os.path.join(db_dir, f"{db_basename}.taxonomy.tsv")

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
        expand(os.path.join(out_dir, "{searchtype}", f"{basename}-x-{db_basename}.{{searchtype}}.{{end}}.tsv"), searchtype = ['blastn', 'diamond-blastx'], end = ['best', 'all']),
        expand(os.path.join(out_dir, f"{basename}-x-{db_basename}.{{end}}.tsv"), end = ['bestblast', 'merged.missed', 'blastn.missed']),


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
        spillover_csv = spillover_csv, 
        tax = tax_file,
    output:
        all = os.path.join(out_dir, "{searchtype}", "{basename}-x-{db_basename}.{searchtype}.all.tsv"),
        best = os.path.join(out_dir, "{searchtype}", "{basename}-x-{db_basename}.{searchtype}.best.tsv"),
    params:
        out_base = os.path.join(out_dir, "{searchtype}", "{basename}-x-{db_basename}.{searchtype}"),
    log: os.path.join(logdir, "select_besthits", "{basename}-x-{db_basename}.{searchtype}.log")
    benchmark: os.path.join(logdir, "select_besthits", "{basename}-x-{db_basename}.{searchtype}.benchmark")
    conda: "conf/env/reports.yml"
    shell:
        """
        python -Werror select-besthits.py \
               --blast {input.blast} \
               --tax {input.tax} \
               --spillover-csv {input.spillover_csv} \
               --output-base {params.out_base} 2> {log}
        """

rule combine_best:
    input:
        blastn = os.path.join(out_dir, "blastn", "{basename}-x-{db_basename}.blastn.best.tsv"),
        blastx = os.path.join(out_dir, "diamond-blastx", "{basename}-x-{db_basename}.diamond-blastx.best.tsv"),
        spillover_csv = spillover_csv,
    output:
        besthit_classif = os.path.join(out_dir, "{basename}-x-{db_basename}.bestblast.tsv"),
        merged_missed = os.path.join(out_dir, "{basename}-x-{db_basename}.merged.missed.tsv"),
        blastn_missed = os.path.join(out_dir, "{basename}-x-{db_basename}.blastn.missed.tsv"),
    params:
        out_base = os.path.join(out_dir, "{basename}-x-{db_basename}"),
    log: os.path.join(logdir, "combine_best", "{basename}-x-{db_basename}.log")
    benchmark: os.path.join(logdir, "combine_best", "{basename}-x-{db_basename}.benchmark")
    conda: 'conf/env/reports.yml'
    shell:
        """
        python -Werror combine-blast.py \
               --spillover-csv {input.spillover_csv} \
               --blastn {input.blastn} \
               --blastx {input.blastx} \
               --output-base {params.out_base} 2> {log}
        """
