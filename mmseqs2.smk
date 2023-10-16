import pandas as pd
import csv
import numpy as np

# download MMSeqs2 with ICTV taxonomy https://zenodo.org/records/6574914/files/virus_tax_db.tar.zst?download=1

basename="2023-08-07-spillover"

onstart:
    print("------------------------------")
    print("mmseqs2 workflow")
    print("------------------------------")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")

rule all:
    input:
        "mmseqs2/virus_tax_db.tar.zst",
        expand(f"output.mmseqs2/{basename}.{{moltype}}", moltype=["dna"]), # "protein"


rule download_mmseqs2_ictv:
    output: "mmseqs2/virus_tax_db.tar.zst"
    shell:
        """
        mkdir -p mmseqs2
        cd mmseqs2
        curl -JLO https://zenodo.org/records/6574914/files/virus_tax_db.tar.zst?download=1
        """

rule mmseqs2_taxonomy_using_dna:
    input: 
        query_fasta = "output.spillover-blast/combined/{basename}.{moltype}.fa.gz",
        database = "mmseqs2/virus_tax_db.tar.zst",
    output: "output.mmseqs2/{basename}.{moltype}"
    conda: "conf/env/mmseqs2.yml"
    shell:
        """
        mmseqs easy-taxonomy {input.query_fasta} {input.ictv_database} {output} tmp --orf-filter 0
        """

# rule mmseqs2_taxonomy_to_tsv:
#     input: 
#         query_fasta = "output.spillover-blast/combined/{basename}.{moltype}.fa.gz",
#         database = "mmseqs2/virus_tax_db.tar.zst",
#     output: "output.mmseqs2/{basename}"
#     conda: "conf/env/mmseqs2.yml"
#     shell:
#         """
#         mmseqs createtsv {input.query_fasta} taxonomyResult taxonomyResult.tsv
#         mmseqs easy-taxonomy {input.query_fasta} {input.ictv_database} {output} tmp --orf-filter 0
#         """

