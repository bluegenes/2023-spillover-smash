import pandas as pd
import csv
import numpy as np

out_dir = "output.sp"
logs_dir = os.path.join(out_dir, 'logs')
basename="spillover"

#sp_file = 'inputs/2023-03-27_spillover_accession-numers.csv'
sp_file = 'inputs/2023-03-27_spillover_accession-numers.head100.csv'

sp = pd.read_csv(sp_file)

null_list = [np.nan, "Unknown", ""]
ACCESSIONS = [a for a in sp['AccessionNumber'] if a and a not in null_list] # don't keep "" entries
# Get a list of non-null and non-"Unknown" values in column1
#non_null_values = vmr.loc[vmr.notnull(vmr["AccessionNumber"]) & (vmr["AccessionNumber"] != "Unknown"), "AccessionNumber"].tolist()
#import pdb;pdb.set_trace()
# print(ACCESSIONS)
#ACCESSIONS = non_null_values

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

rule all:
    input:
        expand(os.path.join(out_dir, f"{basename}.{{moltype}}.zip"), moltype = ['dna', 'protein']),

rule download_spillover_accession:
    output: 
        nucl=protected(os.path.join(out_dir"{basename}/genomic/{acc}.fna.gz")),
        prot=protected(os.path.join("{basename}/protein/{acc}.faa.gz")),
        fileinfo=protected(os.path.join("{basename}/{acc}.fileinfo.csv")),
    conda: "conf/env/biopython.yml"
    log: os.path.join(logs_dir, "downloads", "{basename}/{acc}.log")
    benchmark: os.path.join(logs_dir, "downloads", "{basename}/{acc}.benchmark")
    threads: 1
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="low2",
    shell:
        """
        python genbank_nuccore.py {wildcards.acc} --nucleotide {output.nucl} --protein {output.prot} --fileinfo {output.fileinfo} 2> {log}
        """

rule aggregate_fileinfo_to_fromfile:
    input: 
        fileinfo=expand("{basename}/{acc}.fileinfo.csv", acc=ACCESSIONS, basename="spillover")
    output:
        csv = protected(os.path.join(out_dir, "{basename}.fromfile.csv"))
    run:
        with open(str(output.csv), "w") as outF:
            header = 'name,genome_filename,protein_filename'
            outF.write(','.join(header) + '\n')
            for inp in input:
                with open(str(inp)) as inF:
                    outfile.write(inF.read())

 # Define the checkpoint function that allows us to read the fromfile.csv
checkpoint check_fromfile:
    input: os.path.join(out_dir, f"{basename}.fromfile.csv"),
    output: touch(os.path.join(out_dir,".check_fromfile"))

paramD = {"dna": "dna,k=21,scaled=1,abund", "protein": "protein,k=10,scaled=1,abund"}
rule sketch_fromfile:
    input: 
        fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
        fastas=ancient(Checkpoint_MakePattern("{fn}")),
    output: os.path.join(out_dir, "{basename}.{moltype}.zip")
    params:
        lambda w: paramD[w.moltype]
    threads: 1
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="low2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.log")
    benchmark:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.benchmark")
    shell:
        """
        sourmash sketch fromfile {input.fromfile} -p dna,k=21,scaled=1,abund -o {output} 2> {log}
        """
