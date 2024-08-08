import csv
import pandas as pd
import numpy as np
import urllib

out_dir = "output.vmr"
logs_dir = os.path.join(out_dir, 'logs')

# full VMR
basename = "vmr_MSL39_v1"
vmr_file = 'inputs/VMR_MSL39_v1.acc.bak.tsv'
vmr = pd.read_csv(vmr_file, sep='\t')

# suppressed_records = ['GCF_002987915.1', 'GCF_002830945.1', 'GCF_002828705.1', 'GCA_004789135.1']
# set GenBank Failures dtype to str
# vmr['GenBank Failures'] = vmr['GenBank Failures'].astype(str)
# add 'suppressed' into GenBank Failure column for these records
# vmr.loc[vmr['GenBank Assembly ID'].isin(suppressed_records), 'GenBank Failures'] = 'suppressed'

# null_list = ["", np.nan] + suppressed_records
# ACCESSIONS = [a for a in vmr['GenBank Assembly ID'] if a and a not in null_list] # don't keep "" entries

########################################################################
#### Handle VMR entries with no GenBank Assembly ID or other issues ####
# vmr38_v1: retrieval failures happened for the second of two accs, meaning we still got the valid assembly accession.
# parentheses --> numbers are not the base pair ranges we actually need?
# get working for rest first, then handle parens later
# VMR_ACCESSIONS = []
# vmr_to_genbank = {}
# curate_info = ['suppressed', 'no_assembly', 'multiple_acc'] #, 'parentheses']#, 'retrieval']
# # now select all rows in vmr that had a GenBank Failure
# curate_vmr = vmr.copy()[vmr['GenBank Failures'].isin(curate_info)]

# get clean list of genbank nucleotide accession(s)
def clean_genbank_accs(row):
    genbank_acc = row['Virus GENBANK accession'].split(';')
    # check if ':' is in any and split on that
    for i, acc in enumerate(genbank_acc):
        if ':' in acc:
            genbank_acc[i] = acc.split(':')[1]
        if '(' in acc:
            genbank_acc[i] = acc.split('(')[0]
        genbank_acc[i] = genbank_acc[i].strip()
    return genbank_acc

# if we have rows that need to be curated, find genbank accessions and build vmr_to_genbank dictionary
# if not curate_vmr.empty:
#     curate_vmr['VMR_Accession'] = 'VMR_MSL38_' + curate_vmr['Sort'].astype(str)
#     VMR_ACCESSIONS = curate_vmr['VMR_Accession'].tolist()
#     parentheses_acc = curate_vmr[curate_vmr['GenBank Failures'] == 'parentheses']['VMR_Accession'].tolist()
    
#     curate_vmr.loc[:, 'genbank_accessions'] = curate_vmr.apply(clean_genbank_accs, axis=1)
#     # get dictionary of vmr accession to genbank accession(s)
#     vmr_to_genbank = dict(zip(curate_vmr['VMR_Accession'], curate_vmr['genbank_accessions']))
########################################################################

wildcard_constraints:
    acc = '[^/]+',
    vmr_acc = '[^/]+',

# class Checkpoint_MakePattern:
#     def __init__(self, pattern):
#         self.pattern = pattern
    
#     def get_filenames(self, basename=None, moltype=None):
#         df = pd.read_csv(f"{out_dir}/{basename}.fromfile.csv")
#         filename_col = 'genome_filename'
#         if moltype == "protein":
#             filename_col = 'protein_filename'
#         # filter df to get non-empties in relevant column
#         fastas = df[filename_col][df[filename_col].notnull()].tolist()

#         return fastas

#     def __call__(self, w):
#         global checkpoints
#         # wait for the results of 'check_fromfile'; this will trigger an
#         # exception until that rule has been run.
#         checkpoints.check_fromfile.get(**w)
#         #expand the pattern
#         fastas = self.get_filenames(**w)

#         pattern = expand(self.pattern, fn =fastas,**w)
#         return pattern

rule all:
    input:
        os.path.join(out_dir, f"{basename}.directsketch.zip"),
        os.path.join(out_dir, f"{basename}.nuccore-directsketch.zip"),
        # expand(os.path.join(out_dir, f"{basename}.{{moltype}}.zip"), moltype = ['dna', 'protein']),
        # os.path.join(out_dir, f'{basename}.taxonomy.csv'),
        # os.path.join(out_dir, f'{basename}.protein-taxonomy.csv'),
        # expand(os.path.join(out_dir, f"{basename}.{{moltype}}.lengths.csv"), moltype = ['dna', 'protein']),


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
        ds_csv = os.path.join(out_dir, f"{basename}.directsketch.csv"),
        curated_ds = os.path.join(out_dir, f"{basename}.curate-directsketch.csv"),
        curate_info = os.path.join(out_dir, f"{basename}.ncfasta-to-curate.csv"),
        suppressed = os.path.join(out_dir, f"{basename}.suppressed.csv"),
        lengths = os.path.join(out_dir, f"{basename}.lengths.csv"),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=90,
        time=60,
        partition="low2",
    run:
       # open the assembly summary files and get the "good" (current) and "bad" (historical, suppressed) accessions
        acc2info = {}
        with open(input.good_acc, 'rt') as gacc:
            # good_accs = set()
            r = csv.reader(gacc, delimiter='\t')
            for row in r: 
                if row[0].startswith('#'):
                    continue
                # good_accs.add(row[0])
                acc = row[0]
                gcf_acc = row[0].replace("GCA", "GCF")
                len = row[25] # get genome length! (so we don't need to dl fasta to count :) 
                ftp_path = row[19]
                acc2info[acc] =  (len,ftp_path)
                acc2info[gcf_acc] = (len,ftp_path)
        with open(input.bad_acc, 'rt') as bacc:
            bad_accs = set()
            r = csv.reader(bacc, delimiter='\t')
            for row in r:
                if row[0].startswith('#'):
                    continue
                acc = row[0]
                gcf_acc = row[0].replace("GCA", "GCF")
                bad_accs.add(acc)
                bad_accs.add(gcf_acc)
        # open the vmr file and write directsketch (gbsketch) file for "good" (current) accessions
        curate_ds = open(output.curated_ds, 'w')
        curate_ds.write("accession,name,moltype,md5sum,download_filename,url\n")
        # write a file with info for curation
        curate_info = open(output.curate_info, 'w')
        curate_info.write("name,assembly_failure,genbank_accessions\n")
        lengths = open(output.lengths, 'w')
        lengths.write("accession,length\n")
        with open(vmr_file, 'rt') as infp, open(output.ds_csv, 'w') as outF, open(output.suppressed, 'w') as supF:
            outF.write("accession,name,ftp_path\n")
            r = csv.DictReader(infp, delimiter='\t')
            # write row header to supF
            supF.write(','.join(r.fieldnames) + '\n')
            for row in r:
                fail_reason = None
                acc = row['GenBank Assembly ID']
                if acc and acc not in bad_accs:# and acc in acc2info.keys(): 
                    name = row['GenBank Assembly ID'] + ' ' +  row['Virus name(s)'] + ' ' + row['Virus isolate designation']
                    name = name.strip()
                    if ',' in name:
                        name = name.replace(',', ';')
                    len, ftp_path = acc2info.get(acc, ("", ""))
                    outF.write(f"{acc},{name},{ftp_path}\n")
                    lengths.write(f"{acc},{len}\n")
                else:
                    # this should include ALL failures, including failure to find assembly, suppressed records, etc
                    # write file with "bad" historial accessions so we can curate from their non-assembly dataset records
                    # may not always be desirable, b/c they were likely suppressed for a reason.
                    virus_name = row['Virus name(s)'] + ' ' + row['Virus isolate designation'] #+ '; genbank accs: ' + row["Virus GENBANK accession"]
                    virus_name = virus_name.strip()
                    if ',' in virus_name:
                        virus_name = virus_name.replace(',', ';')
                    gb_acc = row["Virus GENBANK accession"]
                    if acc:
                        print(f"Skipping {acc} for {virus_name} because it is in the historical assembly summary (suppressed)")
                        row['GenBank Failures'] = "suppressed"
                        supF.write(','.join(row.values()) + '\n')
                    else:
                        print(f"Failed to find GenBank Assembly ID for {virus_name}")
                    
                    print(f"using GenBank accessions instead (non-assembly): {gb_acc}")
                    # write separate directsketch (urlsketch) w/ links for genbank accessions that have no assembly
                    for gba in gb_acc.split(';'):
                        gba = gba.strip()
                        if gba:
                            this_name = gba + ' ' + virus_name
                            moltype = "DNA"
                            md5sum = ""
                            dl_filename = f"genbank/nuccore/{gba}.fna.gz"
                            dl_link = f"https://www.ncbi.nlm.nih.gov/nuccore/{gba}?report=fasta"
                            curate_ds.write(f"{gba},{this_name},{moltype},{md5sum},{dl_filename},{dl_link}\n")
                    
                            # create new name for record so we can combine segments if needed
                            vmr_acc = basename + '_'+ str(row['Sort'])
                            curated_name = vmr_acc + ' ' + virus_name
                            curated_name = curated_name.strip()
                            fail_reason = row['GenBank Failures']
                            curate_info.write(f"{curated_name},{fail_reason},{row['Virus GENBANK accession']}\n")

                        # parentheses_acc = curate_vmr[curate_vmr['GenBank Failures'] == 'parentheses']['VMR_Accession'].tolist()
                        # curate_vmr.loc[:, 'genbank_accessions'] = curate_vmr.apply(clean_genbank_accs, axis=1)
                        # get dictionary of vmr accession to genbank accession(s)
                        # vmr_to_genbank = dict(zip(curate_vmr['VMR_Accession'], curate_vmr['genbank_accessions']))
        curate_ds.close()
        curate_info.close()
        lengths.close()

rule directsketch_assembly_datasets:
    input:
        csvfile = os.path.join(out_dir, f"{basename}.directsketch.csv"),
    output:
        zipf = os.path.join(out_dir, f"{basename}.directsketch.zip"),
        failed = os.path.join(out_dir, f"{basename}.dna-directsketch-failed.csv"),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    conda: "conf/env/directsketch.yml"
    log: os.path.join(logs_dir, "directsketch", f"{basename}.log")
    benchmark: os.path.join(logs_dir, "directsketch", f"{basename}.benchmark")
    shell:
        """
        sourmash scripts gbsketch -o {output.zipf} {input.csvfile} \
                                  -p dna,k=21,k=31,scaled=1,abund \
                                  -p protein,k=7,k=10,scaled=1,abund \
                                  --failed {output.failed} 2> {log}
        """
### still need to get FASTA length info -- can I get it from ncbi somehow?
### for curated accessions, we need to download the fasta files from GenBank + find the protein accessions + their FASTA

# download nuccore fasta for curated fasta
rule directsketch_curated:
    input:
        csvfile = os.path.join(out_dir, f"{basename}.curate-directsketch.csv"),
    output:
        zipf = os.path.join(out_dir, f"{basename}.nuccore-directsketch.zip"),
        failed = os.path.join(out_dir, f"{basename}.nuccore-directsketch-failed.csv"),
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    conda: "conf/env/directsketch.yml"
    params:
        fastadir= "genbank/nuccore",
    log: os.path.join(logs_dir, "directsketch", f"{basename}.curate.log")
    benchmark: os.path.join(logs_dir, "directsketch", f"{basename}.curate.benchmark")
    shell:
        """
        sourmash scripts urlsketch -o {output.zipf} {input.csvfile} \
                                  -p dna,k=21,k=31,scaled=1,abund \
                                  -k --fasta {params.fastadir} \
                                  --failed {output.failed} 2> {log}
        """

# rule sigcat_curated:
#     input:
#         sigs = os.path.join(out_dir, f"{basename}.nuccore-directsketch.zip"),
#     output:
#         sig = os.path.join(out_dir, f"{basename}.nuccore-directsketch.sig"),
#     conda: "conf/env/sourmash.yml"
#     shell:
#         """
#         sourmash sig cat {input.sigs} -o {output}
#         """


# # Where we don't have assembly datasets, curate fasta files from GenBank
# rule download_genbank_accession:
#     output: 
#         nucl=protected(os.path.join(out_dir, "genbank/nuccore/nucleotide/{acc}.fna.gz")),
#         fileinfo=protected(os.path.join(out_dir, "genbank/curated/fileinfo/{acc}.fileinfo.csv")),
#     params:
#         prot_dir= os.path.join(out_dir, "genbank/curated/protein"),
#         prot=protected(os.path.join(out_dir, "genbank/curated/protein/{acc}.faa.gz")),
#     conda: "conf/env/biopython.yml"
#     log: os.path.join(logs_dir, "downloads", "{acc}.log")
#     benchmark: os.path.join(logs_dir, "downloads", "{acc}.benchmark")
#     threads: 1
#     resources:
#         mem_mb=3000,
#         runtime=60,
#         time=90,
#         partition="low2",
#     shell:
#         """
#         mkdir -p {params.prot_dir}
#         python genbank_nuccore.py {wildcards.acc} --nucleotide {output.nucl} --protein {params.prot} --fileinfo {output.fileinfo} 2> {log}
#         """

# cat the individual fasta files (genbank accessions) into a single fasta file per vmr accession
# rule cat_components_to_vmr_assembly:
#     input: 
#         fileinfo=lambda w: expand(os.path.join(out_dir, "genbank/curated/fileinfo/{acc}.fileinfo.csv"), acc=vmr_to_genbank[w.vmr_acc])
#     output:
#         nucl= protected(os.path.join(out_dir, "curated/nucleotide/{vmr_acc}.fna.gz")),
#         fileinfo= protected(os.path.join(out_dir, "curated/fileinfo/{vmr_acc}.fileinfo.csv")),
#     params:
#         prot_dir= os.path.join(out_dir, "curated/protein"),
#         prot_out=protected(os.path.join(out_dir, "curated/protein/{vmr_acc}.faa.gz")),
#     log: os.path.join(logs_dir, "curate_fasta", "{vmr_acc}.log")
#     benchmark: os.path.join(logs_dir, "curate_fasta", "{vmr_acc}.benchmark")
#     threads: 1
#     resources:
#         mem_mb=3000,
#         time=90,
#         partition="med2",
#     shell:
#         """
#         mkdir -p {params.prot_dir}
#         python curate-fasta.py --input-fileinfo {input.fileinfo} --curated-acc {wildcards.vmr_acc} \
#                                --nucl-out {output.nucl} --prot-out {params.prot_out} \
#                                --fileinfo-out {output.fileinfo} 2> {log}
#         """

# if genome download failed, download genome fasta so we can translate it
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
