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
curate_info = ['suppressed', 'no_assembly', 'multiple_acc', 'parentheses']#, 'retrieval']
# now select all rows in vmr that had a GenBank Failure
rows_to_curate = vmr[vmr['GenBank Failures'].isin(curate_info)]

# for these, try downloading assembly again
#assembly_accs = rows_to_curate[rows_to_curate['GenBank Failures'] == 'retrieval']['GenBank Assembly ID'].tolist()
#download_assembly = rows_to_curate[rows_to_curate['GenBank Failures'].isin(['retrieval', 'parentheses'])]['GenBank Assembly ID'].tolist()

# for these cases, download the genbank accession, which should still exist. May be multiple? Cat together if need be.
#parentheses = rows_to_curate[rows_to_curate['GenBank Failures'] == 'parentheses']
#multiple_acc = rows_to_curate[rows_to_curate['GenBank Failures'] == 'multiple_acc']
curate_vmr = vmr.copy()[vmr['GenBank Failures'].isin(['no_assembly', 'suppressed', 'multiple_acc', 'parentheses'])]

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


# import pdb;pdb.set_trace()

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
        #os.path.join(out_dir, f"{basename}.curated-fromfile.csv")
        #expand(os.path.join(out_dir, f"{basename}.{{moltype}}.zip"), moltype = ['dna','protein']),
# make non_asd.fromfile.csv, with the VMR identifiers

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
        prot_out=protected(os.path.join(out_dir, "curated/protein/{vmr_acc}.faa.gz")),
        #nucl=lambda w: expand(os.path.join(out_dir, "genbank/curated/nucleotide/{acc}.fna.gz"), acc=vmr_to_genbank[w.vmr_acc]),
        #prot=lambda w: expand(os.path.join(out_dir, "genbank/curated/protein/{acc}.faa.gz"), acc=vmr_to_genbank[w.vmr_acc]),
    log: os.path.join(logs_dir, "zcat_fasta", "{vmr_acc}.log")
    benchmark: os.path.join(logs_dir, "zcat_fasta", "{vmr_acc}.benchmark")
    threads: 1
    resources:
        mem_mb=3000,
        time=90,
        partition="low2",
    run:
        import gzip
        nucl,prot=[],[]
        # read in filenames
        for inFile in input.fileinfo:
            with open(str(inFile)) as inF:
                reader = csv.DictReader(inF)
                for row in reader:
                    nucl.append(row['genome_filename'])
                    prot.append(row['protein_filename'])
        # combine nucl fastas
        with gzip.open(str(output.nucl), 'wt') as outF:
            for nc in nucl:
                # write nc file to outF
                with open(str(nc)) as inF:
                    outF.write(inF.read())
        if prot: #if protein fastas, combine them
            with gzip.open(str(params.prot_out), 'wt') as outF:
                for pt in prot:
                    with open(str(pt)) as inF:
                        outF.write(inF.read())
        # write fileinfo
        with open(str(output.fileinfo), "w") as outF:
            outF.write(f"name,genome_filename,protein_filename\n")
            if prot:
                outF.write(f"{wildcards.vmr_acc},{output.nucl},{params.prot_out}\n")
            else:
                outF.write(f'{wildcards.vmr_acc},{output.nucl},\n')

        # """
        # zcat {input.nucl} > {output.nucl} 2> {log}
        # zcat {params.prot} > {output.prot} 2>> {log}
        # """ 


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
