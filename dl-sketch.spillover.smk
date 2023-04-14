import pandas as pd
import csv

out_dir = "output.vmr"
logs_dir = os.path.join(out_dir, 'logs')

sp_file = 'inputs/2023-03-27_spillover_accession-numers.csv'

vmr = pd.read_csv(sp_file)

#ACCESSIONS = sp['AccessionNumber'].tolist()
#ACCESSIONS = [a for a in ACCESSIONS if a]
ACCESSIONS = [a for a in sp['AccessionNumber'] if a] # don't keep "" entries

rule all:
    input: 
        prot= os.path.join(out_dir, "spillover.protein.zip"),
        dna = os.path.join(out_dir, "spillover.dna.zip"),

rule download_spillover_accession:
    output: 
        nucl="{basename}/genomic/{acc}.fna.gz",
        prot="{basename}/protein/{acc}.faa.gz",
        fileinfo="{basename}/{acc}.fileinfo.csv",
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
        python genbank_nuccore.py {wildcards.accession} --nucleotide {output.nucl} --protein {output.protein} --fileinfo {output.fileinfo} 2> {log}
        """

rule aggregate_fileinfo_to_fromfile:
    input: 
        fileinfo=expand("{basename}/{acc}.fileinfo.csv", acc=ACCESSIONS, basename="spillover")
    output:
        csv = os.path.join(out_dir, "{basename}.fromfile.csv")
    run:
        with open(str(output.csv), "w") as outF:
            header = 'name,genome_filename,protein_filename'
            outF.write(','.join(header) + '\n')
            for inp in input:
                with open(str(inp)) as inF:
                    outfile.write(inF.read())

rule sketch_fromfile_dna:
    input: os.path.join(out_dir, "{basename}.fromfile.csv")
    output: os.path.join(out_dir, "{basename}.dna.zip")
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "sketch", "{basename}.log")
    benchmark:  os.path.join(logs_dir, "sketch", "{basename}.benchmark")
    threads: 1
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="low2",
    shell:
        """
        sourmash sketch fromfile {input} -p dna,k=21,scaled=1,abund -o {output} 2> {log}
        """

rule sketch_fromfile_protein:
    input: os.path.join(out_dir, "{basename}.fromfile.csv")
    output: os.path.join(out_dir, "{basename}.protein.zip")
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "sketch", "{basename}.log")
    benchmark:  os.path.join(logs_dir, "sketch", "{basename}.benchmark")
    threads: 1
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="low2",
    shell:
        """
        sourmash sketch fromfile {input} -p protein,k=10,scaled=1,abund -o {output} 2> {log}
        """