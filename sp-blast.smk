import pandas as pd
import csv

db_dir = "output.vmr"
out_dir = "output.spillover-blast"
logdir = os.path.join(out_dir, 'logs')
basename="sfv"
db_basename = "ictv-spumavirus"

# test simian foamy virus
sm_file = 'inputs/simian-foamy.spillover.csv'
sm = pd.read_csv(sm_file)
ACCS = sm['AccessionNumber']
# get genome filenames
fromfile = 'output.spillover/spillover.fromfile.csv'
ff = pd.read_csv(fromfile)
ff['ident'] = ff['name'].str.split(' ', expand=True)[0]
ff.set_index('ident', inplace=True)

# todo..figure out which protein filenames exist.. (forblastp)

onstart:
    print("------------------------------")
    print("blast workflow")
    print("------------------------------")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")





class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_filenames(self, basename=None, moltype=None):
        df = pd.read_csv(f"{out_dir}/{basename}.fromfile.csv")
        filename_col = 'genome_filename'
        if moltype == "protein":
            filename_col = 'protein_filename'
        # filter df to get non-empties in relevant column
        fastas = df[filename_col][df[filename_col].notnull()].tolist()

        return fastas

    def __call__(self, w):
        global checkpoints
        # wait for the results of 'check_fromfile'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_fromfile.get(**w)
        #expand the pattern
        fastas = self.get_filenames(**w)

        pattern = expand(self.pattern, fn =fastas,**w)
        return pattern

# Define a rule to execute all blast searches
rule all:
    input:
        expand(os.path.join(out_dir, 'blast', "{acc}-x-{db_basename}.blast.txt"), acc = ACCS ,db_basename=db_basename),
        expand(os.path.join(out_dir, 'diamond', "{acc}-x-{db_basename}.diamond-{type}.txt"), acc = ACCS, db_basename=db_basename, type = ['blastx']),#, 'blastp'])

rule blastn:
    input:
        query = lambda w: ff.at[w.acc, 'genome_filename'],
        index = os.path.join(db_dir, "blast", "{db_basename}.dna.index.nhr")
    output:
        results = os.path.join(out_dir, 'blast', "{acc}-x-{db_basename}.blast.txt")
    log: os.path.join(logdir, 'blastn', "{acc}-x-{db_basename}.blastn.log")
    benchmark: os.path.join(logdir, 'blastn', "{acc}-x-{db_basename}.blastn.benchmark")
    conda: "conf/env/blast.yml"
    threads: 4
    shell:
        """
        blastn -query {input.query} -db {input.index} -out {output.results} \ 
               -num_threads {threads} 2> {log}
        """


# Rule to perform DIAMOND blastp
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
                       --out {output.results} --quiet \
                       --sensitive --threads {threads} 2> {log}
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


