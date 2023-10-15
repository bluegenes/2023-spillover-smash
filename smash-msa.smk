import os
import csv
import pandas as pd

# Load the mapping info from the CSV into a DataFrame
# searchDF = pd.read_csv('spillover-x-vmr38.manysearch0.k21.csv')
# gatherDF = pd.read_csv('spillover-x-vmr38.mgather0.with-lineages.csv')

# make dictionary from query_name, match_name columns (aggregate query_name and use as key for all match names)
# sMap =  searchDF.groupby('query_name')['match_name'].apply(list).to_dict()

# for now, let's just load a query csv and a reference csv and map all query --> all ref
basename = 'mammarenavirus'
out_dir = 'output.mm'
query_csv = 'inputs/mm.spillover.csv'
ref_csv = 'inputs/mm.vmr.tsv'
fromfile_csv = f'output.spillover/mammarenavirus.fromfile.csv'
fastaDF = pd.read_csv(fromfile_csv)
spillover_fastas = {'dna': fastaDF['genome_filename'], 'protein': fastaDF['protein_filename']}
rInf = pd.read_csv(ref_csv, sep = '\t')
# basename = 'orthohantavirus'
# out_dir = 'output.orthohantavirus'
# query_csv = 'inputs/hantavirus.spillover.csv'
# ref_csv = 'inputs/orthohantavirus.vmr.csv'
# rInf = pd.read_csv(ref_csv, sep = ',')


qInf = pd.read_csv(query_csv)
QUERY_ACC = qInf['AccessionNumber']
REF_ACC = rInf['GenBank Assembly ID'] # future: need to handle curated / non Assembly datasets identifiers too!
logs_dir = os.path.join(out_dir, 'logs')

rule all:
    input:
        expand(os.path.join(out_dir, "{basename}-x-{ref}.minimap2.bam"), basename = basename, ref = REF_ACC),
        expand(os.path.join(out_dir, "extracted", "{basename}-x-{ref}.extract.fna"), basename=basename, ref = REF_ACC),
        # expand(os.path.join(out_dir, "mafft", "{basename}-x-{ref}.alignment.fasta"), basename=basename, ref= REF_ACC),
        os.path.join(out_dir, "mafft_allref", f"{basename}.alignment.fasta"),
        os.path.join(out_dir, "extracted", f"{basename}.ref-coodinates.csv"),
        os.path.join(out_dir, f"{basename}.classif.csv"),
        expand(os.path.join(out_dir, 'muscle', f'{basename}.{{moltype}}.fa'), moltype = ['protein']),
        # expand("phylogenetic_tree_{query}_{ref}.nwk", query=sMap.keys(), reference=lambda x: mapping_dict[x])

rule combine_dna_fasta:
    input: 
        query_fastas = expand("output.spillover/genomic/{q_acc}.fna.gz", q_acc = QUERY_ACC),
    output: os.path.join(out_dir, "combined", "{basename}.dna.fna")
    log:  os.path.join(logs_dir, "combine", "{basename}.dna.log")
    benchmark:  os.path.join(logs_dir, "combine", "{basename}.dna.benchmark")
    shell:
        """
        gunzip -c {input.fastas} > {output} 2> {log}
        """

rule combine_protein_fasta:
    input: 
        # query_fastas = expand("output.spillover/protein/{q_acc}.faa.gz", q_acc = QUERY_ACC),
        query_fastas = lambda w: spillover_fastas['protein'],
        ref_fastas = expand('genbank/proteomes/{ref}_protein.faa.gz', ref = REF_ACC),
    output: 
        query_combined = os.path.join(out_dir, "combined", "{basename}.protein.faa"),
        query_and_ref = os.path.join(out_dir, "combined", "{basename}.wrefs.protein.faa"),
    log:  os.path.join(logs_dir, "combine", "{basename}.protein.log")
    benchmark:  os.path.join(logs_dir, "combine", "{basename}.protein.benchmark")
    shell:
        """
        gunzip -c {input.query_fastas} > {output.query_combined} 2> {log}
        gunzip -c {input.query_fastas} {input.ref_fastas} > {output.query_and_ref} 2>> {log}
        """

rule samtools_index_reference:
    input: "genbank/genomes/{r_acc}_genomic.fna.gz",
    output:
        fa = "genbank/genomes/{r_acc}_genomic.fna",
        fai = "genbank/genomes/{r_acc}_genomic.fna.fai",
    conda: "minimap2"
    shell:
        """
        gunzip -c {input} > {output.fa}
        samtools faidx {output.fa}
        """ 
# map all queries to all refs??? Or all queries to EACH ref? Start with the latter
rule minimap2_map_sequences: 
    input:
        # query_fna = "output.spillover/genomic/{q_acc}.fna.gz",
        query_fna = os.path.join(out_dir, "combined", "{basename}.dna.fna"),
        ref_fna = "genbank/genomes/{r_acc}_genomic.fna.gz",
        # genbank/genomes/GCF_013087075.1_genomic.fna.gz
    output:
        os.path.join(out_dir, "minimap2", "{basename}-x-{r_acc}.minimap2.bam"),
        # "alignments_{wildcards.query}_{wildcards.reference}.bam"
        # temp("alignments_{wildcards.query}_{wildcards.reference}.bam")
    # conda: 'conf/env/minimap2.yml'
    conda: "minimap2"
    # conda: "/Users/ntward/miniforge3/envs/minimap2"
    shell:
        """
        minimap2 -ax sr {input.ref_fna} {input.query_fna} | samtools view -Sb - > {output}
        """
# rule bwa_map_sequences: # map all queries to all refs??? Or all queries to EACH ref? Start with the latter
#     input:
#         # query_fna = "output.spillover/genomic/{q_acc}.fna.gz",
#         query_fna = os.path.join(out_dir, "combined", "{basename}.dna.fna.gz")
#         ref_fna = "genbank/genomes/{r_acc}.fna.gz",
#     output:
#         os.path.join(out_dir, "{basename}-x-{r_acc}.bwa.bam"),
#     # conda: 'conf/env/bwa.yml'
#     conda: "bwa"
#     shell:
#         """
#         """

rule samtools_sort:
    input: 
        sorted = os.path.join(out_dir, "minimap2","{basename}-x-{r_acc}.minimap2.bam"),
    output: 
        sorted=os.path.join(out_dir, "minimap2", "{basename}-x-{r_acc}.minimap2.sort.bam"),
        bai = os.path.join(out_dir, "minimap2","{basename}-x-{r_acc}.minimap2.sort.bam.bai"),
    log: os.path.join(logs_dir, "samtools_sort", "{basename}-x-{r_acc}.minimap2.sort.log"),
    conda: "minimap2"
    shell:
        """
        samtools sort {input} -o {output.sorted} 2> {log}
        samtools index {output.sorted} 2>> {log}
        """ 

rule extract_genomic_regions:
    input:
        bam=os.path.join(out_dir, "minimap2", "{basename}-x-{r_acc}.minimap2.sort.bam"),
        bai = os.path.join(out_dir, "minimap2","{basename}-x-{r_acc}.minimap2.sort.bam.bai"),
        reference="genbank/genomes/{r_acc}_genomic.fna",
        ref_fai="genbank/genomes/{r_acc}_genomic.fna.fai",
    output:
        fasta=os.path.join(out_dir, "extracted", "{basename}-x-{r_acc}.extract.fna"),
        coords=os.path.join(out_dir, "minimap2", "{basename}-x-{r_acc}.minimap2.consensus_coordinates.txt")
    log: os.path.join(logs_dir, "extract", "{basename}-x-{r_acc}.extract.log")
    shell:
        """
        python find-and-extract-regions.py --bam {input.bam} --fasta {input.reference} \
                                           --output-fasta {output.fasta} --output-coords {output.coords} \
                                           2> {log}
        """

rule align_sequences_single_ref:
    input:
        queries=os.path.join(out_dir, "combined", "{basename}.dna.fna"),
        extracted=os.path.join(out_dir, "extracted", "{basename}-x-{r_acc}.extract.fna"),
    output:
        combined = os.path.join(out_dir, "combined", "{basename}_{r_acc}.fna"),
        alignment=os.path.join(out_dir, "mafft", "{basename}-x-{r_acc}.alignment.fasta"),
    log: os.path.join(logs_dir, "mafft", "{basename}-x-{r_acc}.alignment.log"),
    shell:
        """
        cat {input.queries} {input.extracted} > {output.combined} 2>{log}
        mafft --auto {output.combined} > {output.alignment} 2>> {log}
        """

# ultimately, we need groupings of queries to use with groups of reference sequences. 
# perhaps dictionary of {query_name: closest ref} then group...
# so then [group of queries] : group of refs
# or [group of refs]: group of queries
# ...
rule align_sequences_grouped_refs:
    input:
        queries=os.path.join(out_dir, "combined", "{basename}.dna.fna"),
        extracted=expand(os.path.join(out_dir, "extracted", "{basename}-x-{r_acc}.extract.fna"), basename=basename,r_acc=REF_ACC),
    output:
        combined = os.path.join(out_dir, "combined", "{basename}_w_refs.fna"),
        alignment=os.path.join(out_dir, "mafft_allref", "{basename}.alignment.fasta"),
    log: os.path.join(logs_dir, "mafft", "{basename}.alignment.log"),
    shell:
        """
        cat {input.queries} {input.extracted} > {output.combined}  2> {log}
        mafft --auto {output.combined} > {output.alignment} 2>> {log}
        """

# Aggregate the reference genome coordinates and contig names
rule aggregate_coordinates:
    input:
        coords=expand(os.path.join(out_dir, "minimap2", "{{basename}}-x-{r_acc}.minimap2.consensus_coordinates.txt"), r_acc=REF_ACC),
    output:
        os.path.join(out_dir, "extracted", "{basename}.ref-coodinates.csv")
    run:
        with open(str(output), 'w') as outF:
            writer = csv.writer(outF, delimiter=',')
            for coords_file in input.coords:
                with open(coords_file, 'r') as inF:
                    ref_name = os.path.basename(coords_file).split('-x-')[1].rsplit('.minimap2')[0]
                    reader = csv.reader(inF, delimiter='\t')
                    for row in reader:
                        contig, start, end = row
                        writer.writerow([ref_name, contig, start, end])
                    

rule msa_closest:
    input:
        coords= os.path.join(out_dir, "extracted", "{basename}.ref-coodinates.csv"),
        msa=os.path.join(out_dir, "mafft_allref", "{basename}.alignment.fasta"),
    output:
        classif = os.path.join(out_dir, "{basename}.classif.csv"),
    shell:
        """
        python msa-to-closest-ref.py -m {input.msa} \
                                     -c {input.coords} \
                                     -o {output} \
                                     --top-n 10
        """



rule muscle_protein:
    input:
        fasta = os.path.join(out_dir, "combined", "{basename}.wrefs.protein.faa"), 
        # fastas = lambda w: spillover_fastas[w.moltype],
    output:
        aligned = os.path.join(out_dir, 'muscle', '{basename}.{moltype}.fa')
    log: os.path.join(logs_dir, 'muscle', '{basename}.{moltype}.log')
    shell:
        """
        muscle -align {input} -output {output.aligned} 2> {log}
        """

rule iqtree_muscle:
    input:
    output:
    shell:
        """
        iqtree -s {input} --mem 8G -v -m TEST > {output}
        """


# build VCF from original mapping? ANI to closest ref --> classification. Also get ANI to any other ref it matched
# nucl cutoff for diff species/genus, etc
# tree here is useful as verification --> check branch lengths?

# rule build_tree:
#     input:
#         "aligned_sequences_{wildcards.query}_{wildcards.reference}.fasta"
#     output:
#         "phylogenetic_tree_{wildcards.query}_{wildcards.reference}.nwk"
#     shell:
#         "FastTree {input} > {output}"
