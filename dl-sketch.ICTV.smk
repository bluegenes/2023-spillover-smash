import pandas as pd
import csv

out_dir = "output.vmr"
logs_dir = os.path.join(out_dir, 'logs')

# vmr_file = 'inputs/VMR_21-221122_MSL37.acc.csv'
vmr_file = 'outputs/VMR_21-221122_MSL37.head100.csv'
vmr = pd.read_csv(vmr_file)
basename = "ictv"

ACCESSIONS = vmr['GenBank Assembly ID'].tolist()

rule all:
    input: 
        expand(os.path.join(out_dir, "{basename}.{moltype}.zip"), basename = 'ictv', moltype = ['dna', 'protein']),


### Rules for ICTV GenBank Assemblies:
# download genbank genome details; make an info.csv file for entry.
rule make_genome_info_csv:
    output:
        csvfile = 'genbank/info/{acc}.info.csv'
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=90,
        time=60,
        partition="low2",
    conda: "conf/env/biopython.yml"
    shell:
        """
        python -Werror genbank_genomes.py {wildcards.acc} \
            --output {output.csvfile}
        """

# download actual genomes!
rule download_matching_genome_wc:
    input:
        csvfile = ancient('genbank/info/{acc}.info.csv')
    output:
        genome = "genbank/genomes/{acc}_genomic.fna.gz"
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    # conda: "conf/env/biopython.yml"
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            url = row['genome_url']
            name = row['ncbi_tax_name']

            print(f"downloading genome for acc {acc}/{name} from NCBI...",
                    file=sys.stderr)
            with open(output.genome, 'wb') as outfp:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    outfp.write(content)
                    print(f"...wrote {len(content)} bytes to {output.genome}",
                           file=sys.stderr)


# download proteome (when it exists)
rule download_matching_proteome_wc:
    input:
        csvfile = ancient('genbank/info/{acc}.info.csv')
    output:
        proteome = "genbank/proteomes/{acc}_protein.faa.gz"
    threads: 1
    resources:
         mem_mb=3000,
         runtime=60,
         time=90,
         partition="bml",#"low2", # bml
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            url = row['protein_url']
            name = row['ncbi_tax_name']

            print(f"downloading proteome for acc {acc}/{name} from NCBI...",
                file=sys.stderr)
            with open(output.proteome, 'wb') as outfp:
                #try:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    outfp.write(content)
                    print(f"...wrote {len(content)} bytes to {output.proteome}",
                          file=sys.stderr)
                #except:
                #    shell('touch {output}')


rule build_fromfile_from_assemblyinfo:
    input: 
        info=expand("genbank/info/{acc}.info.csv", acc=ACCESSIONS)
    output:
        csv = os.path.join(out_dir, "{basename}.fromfile.csv")
    run:
        with open(str(output.csv), "w") as outF:
            header = 'name,genome_filename,protein_filename'
            outF.write(','.join(header) + '\n')
            for inp in input:
                with open(str(inp)) as inF:
                    r = csv.DictReader(csv_file)# has fields:  ["acc", "genome_url", "protein_url", "assembly_report_url", "ncbi_tax_name", "ncbi_taxid"] 
                    rows = list(r)
                    assert len(rows) == 1
                    row = rows[0]
                    acc = row['acc']
                    name = acc + " " + row['ncbi_tax_name']
                    gf = "genbank/genomes/{acc}_genomic.fna.gz"
                    pf= ""
                    if row["protein_url"]:
                        # do we want to add prodigal to go nucl --> protein where we don't yet have protein fastas to download?
                        pf = "genbank/proteomes/{acc}_protein.faa.gz"
                    outF.write(f"{name},{gf},{pf}\n")


class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern
    
    def get_filename(w):
        df = pd.read_csv(f"{out_dir}/{w.basename}.fromfile.csv")
        filename_col = 'genome_filename'
        if w.moltype == "protein":
            filename_col = 'protein_filename'
        # filter df to get non-empties in relevant column
        fastas = [p for p in df["protein_filename"] if p] # don't include empty entries
        return fastas

    
    def __call__(self, w):
        global checkpoints
        # wait for the results of 'check_fromfile'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_fromfile.get(**w)
        #expand the pattern
        accs = self.get_acc()

        pattern = expand(self.pattern, **w)
        return pattern
    
 # Define the checkpoint function
checkpoint check_fromfile:
    input: os.path.join(out_dir, f"{basename}.fromfile.csv"),
    output: touch(os.path.join(out_dir,".check_fromfile"))


paramD = {"dna": "dna,k=21,scaled=1,abund", "protein": "protein,k=10,scaled=1,abund"}
rule sketch_fromfile:
    input: 
        #os.path.join(out_dir, "checkpoint.fasta_download"),
        os.path.join(out_dir, "{basename}.fromfile.csv"),
        Checkpoint_MakePattern("{fn}"),
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
        sourmash sketch fromfile {input} -p dna,k=21,scaled=1,abund -o {output} 2> {log}
        """


        # filtered_df = df[df[filename_col] != ""]
        # names = filtered_df["name"].tolist()
        # accs = [n.split(' ')[0] for n in names]
        # return accs
        #protein = [p for p in df["protein_filename"] if p]
        #return {"dna": dna, "protein": protein}