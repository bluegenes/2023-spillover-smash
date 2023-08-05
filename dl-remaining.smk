import csv
import pandas as pd
import numpy as np

out_dir = "output.vmr"
logs_dir = os.path.join(out_dir, 'logs')

basename = "vmr_MSL38_v1"
vmr_file = 'inputs/VMR_MSL38_v1.acc.tsv'
vmr = pd.read_csv(vmr_file, sep='\t')

suppressed_records = ['GCF_002987915.1', 'GCF_002830945.1', 'GCF_002828705.1', 'GCA_004789135.1']
# first, add suppressed info into the VMR

#add 'suppressed' into GenBank Failure column for these records
vmr.loc[vmr['GenBank Assembly ID'].isin(suppressed_records), 'GenBank Failures'] = 'suppressed'

# vmr38_v1: retrieval failures happened for the second of two accs, meaning we still got the valid assembly accession.
# parentheses --> numbers are not the base pair ranges we actually need?
# get working for rest first, then handle parens later
curate_info = ['suppressed', 'no_assembly', 'multiple_acc'] #, 'parentheses']#, 'retrieval']
# now select all rows in vmr that had a GenBank Failure
curate_vmr = vmr.copy()[vmr['GenBank Failures'].isin(curate_info)]

curate_vmr['VMR_Accession'] = 'VMR_MSL38_' + curate_vmr['Sort'].astype(str)
VMR_ACCESSIONS = curate_vmr['VMR_Accession'].tolist()
parentheses_acc = curate_vmr[curate_vmr['GenBank Failures'] == 'parentheses']['VMR_Accession'].tolist()

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

curate_vmr.loc[:, 'genbank_accessions'] = curate_vmr.apply(clean_genbank_accs, axis=1)
# get dictionary of vmr accession to genbank accession(s)
vmr_to_genbank = dict(zip(curate_vmr['VMR_Accession'], curate_vmr['genbank_accessions']))

wildcard_constraints:
    acc = '[^/]+',
    vmr_acc = '[^/]+',

class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern
    
    def get_filenames(self, basename=None, moltype=None):
        df = pd.read_csv(f"{out_dir}/{basename}.curated-fromfile.csv")
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
        checkpoints.check_curated_fromfile.get(**w)
        #expand the pattern
        fastas = self.get_filenames(**w)

        pattern = expand(self.pattern, fn =fastas,**w)
        return pattern

rule all:
    input:
        expand(os.path.join(out_dir, f"{basename}.{{moltype}}.curated.zip"), moltype = ['dna','protein']),
        os.path.join(out_dir, f'{basename}.curated-taxonomy.tsv'),

rule download_genbank_accession:
    output: 
        nucl=protected(os.path.join(out_dir, "genbank/curated/nucleotide/{acc}.fna.gz")),
        fileinfo=protected(os.path.join(out_dir, "genbank/curated/fileinfo/{acc}.fileinfo.csv")),
    params:
        prot_dir= os.path.join(out_dir, "genbank/curated/protein"),
        prot=protected(os.path.join(out_dir, "genbank/curated/protein/{acc}.faa.gz")),
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
        mkdir -p {params.prot_dir}
        python genbank_nuccore.py {wildcards.acc} --nucleotide {output.nucl} --protein {params.prot} --fileinfo {output.fileinfo} 2> {log}
        """

# cat the individual fasta files (genbank accessions) into a single fasta file per vmr accession
rule cat_components_to_vmr_assembly:
    input: 
        fileinfo=lambda w: expand(os.path.join(out_dir, "genbank/curated/fileinfo/{acc}.fileinfo.csv"), acc=vmr_to_genbank[w.vmr_acc])
    output:
        nucl= protected(os.path.join(out_dir, "curated/nucleotide/{vmr_acc}.fna.gz")),
        fileinfo= protected(os.path.join(out_dir, "curated/fileinfo/{vmr_acc}.fileinfo.csv")),
    params:
        prot_dir= os.path.join(out_dir, "curated/protein"),
        prot_out=protected(os.path.join(out_dir, "curated/protein/{vmr_acc}.faa.gz")),
    log: os.path.join(logs_dir, "curate_fasta", "{vmr_acc}.log")
    benchmark: os.path.join(logs_dir, "curate_fasta", "{vmr_acc}.benchmark")
    threads: 1
    resources:
        mem_mb=3000,
        time=90,
        partition="med2",
    shell:
        """
        mkdir -p {params.prot_dir}
        python curate-fasta.py --input-fileinfo {input.fileinfo} --curated-acc {wildcards.vmr_acc} \
                               --nucl-out {output.nucl} --prot-out {params.prot_out} \
                               --fileinfo-out {output.fileinfo} 2> {log}
        """


rule aggregate_fileinfo_to_fromfile:
   input: 
       fileinfo=expand(os.path.join(out_dir, "curated/fileinfo/{acc}.fileinfo.csv"), acc=VMR_ACCESSIONS)
   output:
       csv = protected(os.path.join(out_dir, "{basename}.curated-fromfile.csv"))
   run:
       with open(str(output.csv), "w") as outF:
           header = 'name,genome_filename,protein_filename'
           outF.write(header + '\n')
           for inp in input:
               with open(str(inp)) as inF:
                   outF.write(inF.read())

#  Define the checkpoint function that allows us to read the fromfile.csv
checkpoint check_curated_fromfile:
    input: os.path.join(out_dir, f"{basename}.curated-fromfile.csv"),
    output: touch(os.path.join(out_dir,".check_curated_fromfile"))

paramD = {"dna": "dna,k=21,k=31,scaled=1,abund", "protein": "protein,k=7,k=10,scaled=1,abund"}
rule sketch_fromfile:
    input: 
        fromfile=os.path.join(out_dir, "{basename}.curated-fromfile.csv"),
        fastas=ancient(Checkpoint_MakePattern("{fn}")),
    output: os.path.join(out_dir, "{basename}.{moltype}.curated.zip")
    params:
        lambda w: paramD[w.moltype]
    threads: 1
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="low2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.curated.log")
    benchmark:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.benchmark")
    shell:
        """
        sourmash sketch fromfile {input.fromfile} -p {params} -o {output} --report-duplicated --ignore-missing 2> {log}
        """


rule build_dna_taxonomy:
    input:
        vmr_file = vmr_file,
        fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
        curated_fromfile=os.path.join(out_dir, "{basename}.curated-fromfile.csv"),
        fastas=ancient(Checkpoint_MakePattern("{fn}")),
    output:
        tax = os.path.join(out_dir, '{basename}.curated-taxonomy.tsv'),
#        prot_tax = os.path.join(out_dir, '{basename}.protein-taxonomy.csv'), #columns prot_name, dna_acc, full_lineage
    params:
        suppressed_records = ' '.join(suppressed_records),
    conda: 'conf/env/reports.yml'
    log:  os.path.join(logs_dir, "build_taxonomy", "{basename}.curated.log")
    benchmark:  os.path.join(logs_dir, "build_taxonomy", "{basename}.curated.benchmark")
    shell:
        """
        python -Werror make-tax.py --vmr-tsv {input.vmr_file} \
                                   --fromfile {input.fromfile} {input.curated_fromfile} \
                                   --output {output.tax} \
                                   --suppressed-records {params.suppressed_records} 2> {log}
        """
