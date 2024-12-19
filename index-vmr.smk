import csv
import pandas as pd
import numpy as np
import urllib

out_dir = "output.vmr"
logs_dir = os.path.join(out_dir, 'logs')

# full VMR
basename = "vmr_MSL39_v4"
 # first use find-assembly-accessions.py to generate this file.
vmr_file = 'inputs/VMR_MSL39.v4_20241106.acc.tsv'
vmr = pd.read_csv(vmr_file, sep='\t')

wildcard_constraints:
    acc = '[^/]+',
    vmr_acc = '[^/]+',

rule all:
    input:
        os.path.join(out_dir, f"{basename}.gbsketch.zip"),
        os.path.join(out_dir, f"{basename}.curate-ds.zip"),
        os.path.join(out_dir, f"{basename}.combined.zip"),


### Rules for ICTV GenBank Assemblies:
rule download_assembly_summary:
    output:
        summary = 'genbank/assembly_summary.viral.txt',
        historical = 'genbank/assembly_summary_historical.txt', # includes suppressed accessions, etc
    shell: 
        """
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt > {output.summary}
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary_historical.txt > {output.historical} 
        """

rule acc_to_directsketch:
    input:
        vmr_file = vmr_file,
        good_acc = 'genbank/assembly_summary.viral.txt',
        bad_acc = 'genbank/assembly_summary_historical.txt',
    output:
        ds_csv = os.path.join(out_dir, f"{basename}.gbsketch.csv"),
        curated_ds = os.path.join(out_dir, f"{basename}.curate-urlsketch.csv"),
        curate_info = os.path.join(out_dir, f"{basename}.ncfasta-to-curate.csv"),
        suppressed = os.path.join(out_dir, f"{basename}.suppressed.csv"),
        lengths = os.path.join(out_dir, f"{basename}.lengths.csv"),
    shell:
        """
        python acc-to-directsketch.py \
            --vmr_file {input.vmr_file} \
            --good_acc {input.good_acc} \
            --bad_acc {input.bad_acc} \
            --ds_csv {output.ds_csv} \
            --curated_ds {output.curated_ds} \
            --curate_info {output.curate_info} \
            --suppressed {output.suppressed} \
            --lengths {output.lengths} \
            --basename {basename}
        """

rule directsketch_assembly_datasets:
    input:
        csvfile = os.path.join(out_dir, f"{basename}.gbsketch.csv"),
    output:
        zipf = os.path.join(out_dir, f"{basename}.gbsketch.zip"),
        failed = os.path.join(out_dir, f"{basename}.gbsketch-assemblies-failed.csv"),
        ch_failed = os.path.join(out_dir, f"{basename}.gbsketch-assemblies-checksum-failed.csv"),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    conda: "directsketch"
    #conda: "conf/env/directsketch.yml"
    log: os.path.join(logs_dir, "directsketch", f"{basename}.gbsketch.log")
    benchmark: os.path.join(logs_dir, "directsketch", f"{basename}.gbsketch.benchmark")
    shell:
        """
        sourmash scripts gbsketch -o {output.zipf} {input.csvfile} \
                                  -p dna,k=21,k=31,scaled=5 \
                                  -p skipm2n3,k=15,k=17,k=19,k=21,scaled=5 \
                                  -p skipm2n3,k=15,k=17,k=19,k=21,k=23,k=25,scaled=5 \
                                  --failed {output.failed} 2> {log}
        """
                                  #-p protein,k=7,k=10,scaled=1,abund \

### still need to get FASTA length info -- can I get it from ncbi somehow?
### for curated accessions, we need to download the fasta files from GenBank + find the protein accessions + their FASTA

# # Where we don't have assembly datasets, curate fasta files from GenBank nuccore
# here, we use directsketch to merge the component fastas into a single sketch
rule directsketch_curated:
    input: 
        os.path.join(out_dir, f"{basename}.curate-urlsketch.csv"),
    output:
        zipf = os.path.join(out_dir, f"{basename}.curate-ds.zip"), 
        failed = os.path.join(out_dir, f"{basename}.curate-ds-failed.csv"),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    conda: "directsketch"
    #conda: "conf/env/directsketch.yml"
    log: os.path.join(logs_dir, "ds-curate", f"{basename}.curate.log")
    benchmark: os.path.join(logs_dir, "ds-curate", f"{basename}.curate.benchmark")
    shell:
        """
        sourmash scripts urlsketch -o {output.zipf} {input} \
                                  -p dna,k=21,k=31,scaled=5 \
                                  -p skipm2n3,k=15,k=17,k=19,k=21,scaled=5 \
                                  -p skipm2n3,k=15,k=17,k=19,k=21,k=23,k=25,scaled=5 \
                                  --failed {output.failed} 2> {log}
        """
                                  #-p protein,k=7,k=10,scaled=100 \

rule combine_sigs:
    input:
        directsketch = os.path.join(out_dir, f"{basename}.gbsketch.zip"),
        curated = os.path.join(out_dir, f"{basename}.curate-ds.zip"),
    output:
        combined = os.path.join(out_dir, f"{basename}.combined.zip"),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    conda: "directsketch"
    #conda: "conf/env/branchwater.yml"
    log: os.path.join(logs_dir, "combine-sketches", f"{basename}.log")
    benchmark: os.path.join(logs_dir, "combine-sketches", f"{basename}.benchmark")
    shell:
        """
        sourmash sig cat {input.directsketch} {input.curated} -o {output.combined} 2> {log}
        """

# rule build_taxonomy:
#     input:
#         vmr_file = vmr_file,
#         directsketch=os.path.join(out_dir, f"{basename}.directsketch.csv"),
#         curated=os.path.join(out_dir, f"curated/{basename}.curate.fromfile.csv"),
#         combined=os.path.join(out_dir, f"{basename}.combined.zip"),
#     output:
#         tax = os.path.join(out_dir, f'{basename}.taxonomy.tsv'),
#         prot_tax = os.path.join(out_dir, f'{basename}.protein-taxonomy.csv'), #columns prot_name, dna_acc, full_lineage
#     params:
#         suppressed_records = ' '.join(suppressed_records),
#     conda: 'conf/env/reports.yml'
#     log:  os.path.join(logs_dir, "build_taxonomy", f"{basename}.log")
#     benchmark:  os.path.join(logs_dir, "build_taxonomy", f"{basename}.benchmark")
#     shell:
#         """
#         python -Werror make-tax.py --vmr-tsv {input.vmr_file} \
#                                    --fromfile {input.fromfile} \
#                                    --output {output.tax} \
#                                    --suppressed-records {params.suppressed_records} 2> {log}
#         """

## RULES FOR PROTEINS

# if proteome download failed, download genome fasta so we can translate it
# rule directsketch_get_fasta_for_failures:
#     input:
#         csvfile = os.path.join(out_dir, f"{basename}.directsketch-failed.csv"),
#     output:
#         zipf = os.path.join(out_dir, f"{basename}.zip"),
#         failed = os.path.join(out_dir, f"{basename}.prot-directsketch-failed.csv"),
#     threads: 1
#     resources:
#         mem_mb=3000,
#         disk_mb=5000,
#         runtime=60,
#         time=90,
#         partition="low2",
#     conda: "conf/env/directsketch.yml"
#     shell:
#         """
#         sourmash scripts gbsketch -o {output.zipf} {input.csvfile} \
#                                   -p protein,k=7,k=10,scaled=1,abund \
#                                   --failed {output.failed} 2> {logs_dir}/directsketch.log
#         """ 

# rule translate_proteomes:
#     input:
#         csv = os.path.join(out_dir, f"{basename}.fromfile.csv"),
#         nucl = "genbank/genomes/{acc}_genomic.fna.gz",
#     output:
#         translated=protected("genbank/translate/{acc}_protein.faa.gz"),
#     conda: "conf/env/seqkit.yml"
#     log: os.path.join(logs_dir, 'seqkit-translate', '{acc}.log')
#     threads: 1
#     shell:
#         """
#         seqkit translate {input.nucl} | gzip > {output} 2> {log}
#         """

# rule translate_curated:
#     input:
#         csv = os.path.join(out_dir, f"{basename}.fromfile.csv"),
#         nucl = f"{out_dir}/curated/nucleotide/{{acc}}.fna.gz",
#     output:
#         translated=protected(os.path.join(out_dir, 'curated', "translate/{acc}.faa.gz")),
#     conda: "conf/env/seqkit.yml"
#     log: os.path.join(logs_dir, 'seqkit-translate', 'curated', '{acc}.log')
#     threads: 1
#     shell:
#         """
#         seqkit translate {input.nucl} | gzip > {output} 2> {log}
#         """


# localrules: build_fromfile_from_assemblyinfo
# rule build_fromfile_from_assemblyinfo: # include curated fileinfo
#     input: 
#         info=expand("genbank/info/{acc}.info.csv", acc=ACCESSIONS),
#         curated=expand(os.path.join(out_dir, "curated/fileinfo/{acc}.fileinfo.csv"), acc=VMR_ACCESSIONS),
#     output:
#         csv = os.path.join(out_dir, "{basename}.fromfile.csv")
#     threads: 1
#     run:
#         with open(str(output.csv), "w") as outF:
#             header = ["name","genome_filename","protein_filename"]
#             outF.write(','.join(header) + '\n')
#             for inp in input.info:
#                 with open(str(inp)) as csvfile:
#                     csv_reader = csv.DictReader(csvfile)# .# has fields:  ["acc", "genome_url", "protein_url", "assembly_report_url", "ncbi_tax_name", "ncbi_taxid"]
#                     rows = list(csv_reader)
#                     assert len(rows) == 1
#                     row = rows[0]
#                     acc = row['acc']
#                     name = acc + " " + row['ncbi_tax_name']
#                     gf = f"genbank/genomes/{acc}_genomic.fna.gz"
#                     pf= ""
#                     if row["protein_url"]:
#                         # do we want to add prodigal to go nucl --> protein where we don't yet have protein fastas to download?
#                         pf = f"genbank/proteomes/{acc}_protein.faa.gz"
#                     else:
#                         pf = f"genbank/translate/{acc}_protein.faa.gz"
#                     outF.write(f"{name},{gf},{pf}\n")
#             for inp in input.curated:
#                with open(str(inp)) as inF:
#                     name,dna,prot = inF.read().split(',')
#                     this_acc = name.split(' ')[0]
#                     if this_acc not in prot:
#                         prot = f"{out_dir}/curated/translate/{this_acc}.faa.gz\n"    
#                     outF.write(f"{name},{dna},{prot}")
#                 #    outF.write(inF.read())

    
 # Define the checkpoint function that allows us to read the fromfile.csv
# checkpoint check_fromfile:
#     input: os.path.join(out_dir, f"{basename}.fromfile.csv"),
#     output: touch(os.path.join(out_dir,".check_fromfile"))


# # paramD = {"dna": "dna,k=21,k=31,scaled=1,abund", "protein": "protein,k=7,k=10,scaled=1,abund"}
# paramD = {"dna": "dna,k=9,k=11,k=13,k=15,k=17,k=19,k=21,k=31,scaled=1,abund", "protein": "protein,k=3,k=4,k=5,k=6,k=7,k=8,k=9,k=10,scaled=1,abund"}
# rule sketch_fromfile:
#     input: 
#         fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
#         fastas=ancient(Checkpoint_MakePattern("{fn}")),
#     output: os.path.join(out_dir, "{basename}.{moltype}.zip")
#     params:
#         lambda w: paramD[w.moltype]
#     threads: 1
#     resources:
#         mem_mb=6000,
#         runtime=4000,
#         time=4000,
#         partition="high2",
#     conda: "conf/env/branchwater.yml"
#     #conda: "pyo3-branch"
#     log:  os.path.join(logs_dir, "manysketch", "{basename}.{moltype}.log")
#     benchmark:  os.path.join(logs_dir, "manysketch", "{basename}.{moltype}.benchmark")
#     shell:
#         """
#         sourmash scripts manysketch {input.fromfile} -p {params} \
#                                     -o {output} 2> {log}
#         """
#         # sourmash sketch fromfile {input.fromfile} -p {params} -o {output} --report-duplicated --ignore-missing 2> {log}

# rule combine_fasta:
#     input: 
#         fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
#         fastas=ancient(Checkpoint_MakePattern("{fn}")),
#     output: os.path.join(out_dir, "{basename}.{moltype}.fa.gz")
#     params:
#     log:  os.path.join(logs_dir, "combine", "{basename}.{moltype}.log")
#     benchmark:  os.path.join(logs_dir, "combine", "{basename}.{moltype}.benchmark")
#     shell:
#         """
#         gunzip -c {input.fastas} | gzip > {output} 2> {log}
#         """

# rule get_fasta_lengths:
#     input: os.path.join(out_dir, "{basename}.{moltype}.fa.gz")
#     output: os.path.join(out_dir, "{basename}.{moltype}.lengths.csv")
#     run:
#         import screed
#         with screed.open(input[0]) as seqfile, open(output[0], 'w') as outfile:
#             outfile.write("sequence_name,length\n")
#             for record in seqfile:
#                 bp = len(record.sequence)
#                 outfile.write(f"{record.name},{bp}\n")


# # Rule to build BLAST index for the combined gzipped fasta file
# rule build_nucl_index:
#     input:
#         fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
#         fasta = os.path.join(out_dir, "{basename}.dna.fa.gz")
#     output:
#         index = os.path.join(out_dir, "blast", "{basename}.dna.index.nhr")
#     params:
#         title = os.path.join("{basename}"),
#         out_base = os.path.join(out_dir, "blast", "{basename}.dna.index")
#     log:  os.path.join(logs_dir, "blast-index", "{basename}.dna.log")
#     benchmark:  os.path.join(logs_dir, "blast-index", "{basename}.dna.benchmark")
#     conda: "conf/env/blast.yml"
#     shell:
#         """
#         gunzip -c {input.fasta} | makeblastdb -in - -dbtype nucl -parse_seqids \
#                -out {params.out_base} -title {params.title} 2> {log}
#         """

# rule build_prot_index:
#     input:
#         fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
#         fasta = os.path.join(out_dir, "{basename}.protein.fa.gz"),
#     output:
#         index = os.path.join(out_dir, "diamond", "{basename}.protein.fa.gz" + ".dmnd"),
#     log:  os.path.join(logs_dir, "diamond-index", "{basename}.log")
#     benchmark:  os.path.join(logs_dir, "diamond-index", "{basename}.benchmark")
#     conda: "conf/env/diamond.yml"
#     shell:
#         """
#         diamond makedb --in {input.fasta} --db {output.index} --quiet 2> {log}
#         """


# rule build_dna_taxonomy:
#     input:
#         vmr_file = vmr_file,
#         fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
#         fastas=ancient(Checkpoint_MakePattern("{fn}")),
#     output:
#         tax = os.path.join(out_dir, '{basename}.taxonomy.tsv'),
# #        prot_tax = os.path.join(out_dir, '{basename}.protein-taxonomy.csv'), #columns prot_name, dna_acc, full_lineage
#     params:
#         suppressed_records = ' '.join(suppressed_records),
#     conda: 'conf/env/reports.yml'
#     log:  os.path.join(logs_dir, "build_taxonomy", "{basename}.log")
#     benchmark:  os.path.join(logs_dir, "build_taxonomy", "{basename}.benchmark")
#     shell:
#         """
#         python -Werror make-tax.py --vmr-tsv {input.vmr_file} \
#                                    --fromfile {input.fromfile} \
#                                    --output {output.tax} \
#                                    --suppressed-records {params.suppressed_records} 2> {log}
#         """
