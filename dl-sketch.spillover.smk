import pandas as pd
import csv

out_dir = "output.vmr"
logs_dir = os.path.join(out_dir, 'logs')

sp_file = 'inputs/2023-03-27_spillover_accession-numers.csv'

vmr = pd.read_csv(sp_file)

ACCESSIONS = sp['AccessionNumber'].tolist()

rule all:
    input: 
        prot= os.path.join(out_dir, "spillover.protein.zip"),
        dna = os.path.join(out_dir, "spillover.dna.zip"),

rule download_spillover_accession:
    output: 
        nucl="spillover/genomic/{acc}.fna.gz",
        prot="spillover/protein/{acc}.faa.gz",
        fileinfo="spillover/{acc}.fileinfo.csv",
    log: os.path.join(logs_dir, "downloads", "spillover/{acc}.log")
    benchmark: os.path.join(logs_dir, "downloads", "spillover/{acc}.benchmark")
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
        fileinfo=expand("spillover/{acc}.fileinfo.csv", acc=ACCESSIONS)
    output:
        csv = os.path.join(out_dir, "spillover.fromfile.csv")
    run:
        with open(str(output.csv), "w") as outF:
            header = 'name,genome_filename,protein_filename'
            outF.write(','.join(header) + '\n')
            for inp in input:
                with open(str(inp)) as inF:
                    outfile.write(inF.read())

rule sketch_fromfile_dna:
    input: os.path.join(out_dir, "spillover.fromfile.csv")
    output: os.path.join(out_dir, "spillover.dna.zip")
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "sketch", "spillover.log")
    benchmark:  os.path.join(logs_dir, "sketch", "spillover.benchmark")
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
    input: os.path.join(out_dir, "spillover.fromfile.csv")
    output: os.path.join(out_dir, "spillover.protein.zip")
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "sketch", "spillover.log")
    benchmark:  os.path.join(logs_dir, "sketch", "spillover.benchmark")
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