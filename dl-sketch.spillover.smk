import pandas as pd
import csv
import numpy as np

out_dir = "output.spillover"
logs_dir = os.path.join(out_dir, 'logs')
basename="spillover"

# sp_file = 'inputs/2023-03-27_spillover_accession-numers.csv'
sp_file = 'inputs/2023-08-07_spillover_accession-numbers_info.csv'
# mammarenavirus genus only
# sp_file = 'inputs/mm.spillover.csv'
# basename = 'mammarenavirus'
# orthohantavirus only
#sp_file = 'inputs/hantavirus.spillover.csv'
#basename = 'hantavirus'

sp = pd.read_csv(sp_file)

null_list = [np.nan, "Unknown", ""]
ACCESSIONS = [a for a in sp['AccessionNumber'] if a and a not in null_list] # don't keep "" entries

# to test, let's only use the first 1000
#basename="spillover1000"
#ACCESSIONS = ACCESSIONS[:1000]

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
        expand(os.path.join(out_dir, f"{basename}.{{moltype}}.zip"), moltype = ['dna','protein']),

rule download_spillover_accession:
    output: 
        nucl=protected(os.path.join(out_dir, "genomic/{acc}.fna.gz")),
        fileinfo=protected(os.path.join(out_dir, "fileinfo/{acc}.fileinfo.csv")),
    params:
        prot=protected(os.path.join(out_dir, "protein/{acc}.faa.gz")),
    conda: "conf/env/biopython.yml"
    log: os.path.join(logs_dir, "downloads", "{acc}.log")
    benchmark: os.path.join(logs_dir, "downloads", "{acc}.benchmark")
    threads: 1
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="low2",
    shell:
        """
        python genbank_nuccore.py {wildcards.acc} --nucleotide {output.nucl} --protein {params.prot} --fileinfo {output.fileinfo} 2> {log}
        """

rule translate_protein:
    input:
        fileinfo=os.path.join(out_dir, "fileinfo/{acc}.fileinfo.csv"),
        nucl=os.path.join(out_dir, "genomic/{acc}.fna.gz"),
    output:
        translated=protected(os.path.join(out_dir, "translate/{acc}.faa.gz")),
    conda: "conf/env/seqkit.yml"
    log: os.path.join(logs_dir, 'seqkit', '{acc}.log')
    threads: 1
    shell:
        """
        seqkit translate {input.nucl} | gzip > {output} 2> {log}
        """

rule aggregate_fileinfo_to_fromfile:
    input: 
        fileinfo=expand(os.path.join(out_dir, "fileinfo/{acc}.fileinfo.csv"), acc=ACCESSIONS)
    output:
        csv = protected(os.path.join(out_dir, "{basename}.fromfile.csv"))
    run:
        with open(str(output.csv), "w") as outF:
            header = 'name,genome_filename,protein_filename'
            outF.write(header + '\n')
            for inp in input:
                with open(str(inp)) as inF:
                    # outF.write(inF.read())
                    # if protein_filename doesn't exist, use translated
                    name,dna,prot = inF.read().split(',')
                    this_acc = name.split(' ')[0]
                    if this_acc not in prot:
                        prot = f"{out_dir}/translate/{this_acc}.faa.gz\n"
                    outF.write(f"{name},{dna},{prot}")

# Define the checkpoint function that allows us to read the fromfile.csv
checkpoint check_fromfile:
    input: os.path.join(out_dir, f"{basename}.fromfile.csv"),
    output: touch(os.path.join(out_dir,".check_fromfile"))

paramD = {"dna": "dna,k=9,k=11,k=13,k=15,k=17,k=19,k=21,k=31,scaled=1,abund", "protein": "protein,k=5,k=6,k=7,k=8,k=9,k=10,scaled=1,abund"}
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
    # conda: "conf/env/sourmash.yml"
    conda: "conf/env/branchwater.yml"
    log:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.log")
    benchmark:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.benchmark")
    shell:
        """
        sourmash scripts manysketch {input.fromfile} -p {params} \
                                    -o {output} 2> {log}
        """
        # sourmash sketch fromfile {input.fromfile} -p {params} -o {output} --report-duplicated --ignore-missing 2> {log}
