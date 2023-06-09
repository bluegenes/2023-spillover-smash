import csv
import pandas as pd
import numpy as np

out_dir = "output.vmr"
logs_dir = os.path.join(out_dir, 'logs')


basename = "vmr_MSL38_v1"
#vmr_file = 'inputs/VMR_MSL38_v1.acc.csv'
#vmr = pd.read_csv(vmr_file)
vmr_file = 'inputs/VMR_MSL38_v1.acc.tsv'
vmr = pd.read_csv(vmr_file, sep='\t')

suppressed_records = ['GCF_002987915.1', 'GCF_002830945.1', 'GCF_002828705.1', 'GCA_004789135.1']

# subsets and tests
#vmr_file = 'inputs/VMR_MSL38_v1.lassa.tsv'
#vmr = pd.read_csv(vmr_file, sep='\t')
#basename = "vmr_MSL38_v1.lassa"
#vmr_file = 'inputs/spumavirus.VMR_MSL38v1.acc.csv'
#vmr_file = 'inputs/VMR_MSL38_v1.acc.head100.csv'
#vmr_file = 'inputs/VMR_21-221122_MSL37.acc.csv'
#vmr_file = 'outputs/VMR_21-221122_MSL37.head100.csv'
#basename = "ictv-h100"
#basename = "ictv-spumavirus"

null_list = ["", np.nan] + suppressed_records
ACCESSIONS = [a for a in vmr['GenBank Assembly ID'] if a and a not in null_list] # don't keep "" entries
# print(ACCESSIONS)


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
        expand(os.path.join(out_dir, f"{basename}.{{moltype}}.zip"), moltype = ['dna', 'protein']),
        expand(os.path.join(out_dir, "blast", f"{basename}.dna.index.nhr")),
        expand(os.path.join(out_dir, "diamond", f"{basename}.protein.fa.gz.dmnd")),
        os.path.join(out_dir, f'{basename}.taxonomy.csv'),
        os.path.join(out_dir, f'{basename}.protein-taxonomy.csv'),

### Rules for ICTV GenBank Assemblies:
# download genbank genome details; make an info.csv file for entry.
rule make_genome_info_csv:
    output:
        csvfile = 'genbank/info/{acc}.info.csv',
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
        genome = protected("genbank/genomes/{acc}_genomic.fna.gz")
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
        proteome = protected("genbank/proteomes/{acc}_protein.faa.gz")
    threads: 1
    resources:
         mem_mb=3000,
         runtime=60,
         time=90,
         partition="low2", # bml
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

localrules: build_fromfile_from_assemblyinfo
rule build_fromfile_from_assemblyinfo:
    input: 
        info=expand("genbank/info/{acc}.info.csv", acc=ACCESSIONS)
    output:
        csv = os.path.join(out_dir, "{basename}.fromfile.csv")
    threads: 1
    run:
        with open(str(output.csv), "w") as outF:
            header = ["name","genome_filename","protein_filename"]
            outF.write(','.join(header) + '\n')
            for inp in input:
                with open(str(inp)) as csvfile:
                    csv_reader = csv.DictReader(csvfile)# .# has fields:  ["acc", "genome_url", "protein_url", "assembly_report_url", "ncbi_tax_name", "ncbi_taxid"]
                    rows = list(csv_reader)
                    assert len(rows) == 1
                    row = rows[0]
                    acc = row['acc']
                    name = acc + " " + row['ncbi_tax_name']
                    gf = f"genbank/genomes/{acc}_genomic.fna.gz"
                    pf= ""
                    if row["protein_url"]:
                        # do we want to add prodigal to go nucl --> protein where we don't yet have protein fastas to download?
                        pf = f"genbank/proteomes/{acc}_protein.faa.gz"
                    outF.write(f"{name},{gf},{pf}\n")

    
 # Define the checkpoint function that allows us to read the fromfile.csv
checkpoint check_fromfile:
    input: os.path.join(out_dir, f"{basename}.fromfile.csv"),
    output: touch(os.path.join(out_dir,".check_fromfile"))


paramD = {"dna": "dna,k=21,k=31,scaled=1,abund", "protein": "protein,k=7,k=10,scaled=1,abund"}
rule sketch_fromfile:
    input: 
        fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
        fastas=ancient(Checkpoint_MakePattern("{fn}")),
    output: os.path.join(out_dir, "{basename}.{moltype}.zip")
    params:
        lambda w: paramD[w.moltype]
    threads: 1
    resources:
        mem_mb=6000,
        runtime=4000,
        time=4000,
        partition="high2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.log")
    benchmark:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.benchmark")
    shell:
        """
        sourmash sketch fromfile {input.fromfile} -p {params} -o {output} --report-duplicated --ignore-missing 2> {log}
        """

rule combine_fasta:
    input: 
        fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
        fastas=ancient(Checkpoint_MakePattern("{fn}")),
    output: os.path.join(out_dir, "{basename}.{moltype}.fa.gz")
    params:
    log:  os.path.join(logs_dir, "combine", "{basename}.{moltype}.log")
    benchmark:  os.path.join(logs_dir, "combine", "{basename}.{moltype}.benchmark")
    shell:
        """
        zcat {input.fastas} | gzip > {output} 2> {log}
        """


# Rule to build BLAST index for the combined gzipped fasta file
rule build_nucl_index:
    input:
        fasta = os.path.join(out_dir, "{basename}.dna.fa.gz")
    output:
        index = os.path.join(out_dir, "blast", "{basename}.dna.index.nhr")
    params:
        title = os.path.join("{basename}"),
        out_base = os.path.join(out_dir, "blast", "{basename}.dna.index")
    log:  os.path.join(logs_dir, "blast-index", "{basename}.dna.log")
    benchmark:  os.path.join(logs_dir, "blast-index", "{basename}.dna.benchmark")
    conda: "conf/env/blast.yml"
    shell:
        """
        gunzip -c {input.fasta} | makeblastdb -in - -dbtype nucl -parse_seqids \
               -out {params.out_base} -title {params.title} 2> {log}
        """

rule build_prot_index:
    input:
        fasta = os.path.join(out_dir, "{basename}.protein.fa.gz"),
    output:
        index = os.path.join(out_dir, "diamond", "{basename}.protein.fa.gz" + ".dmnd"),
    log:  os.path.join(logs_dir, "diamond-index", "{basename}.log")
    benchmark:  os.path.join(logs_dir, "diamond-index", "{basename}.benchmark")
    conda: "conf/env/diamond.yml"
    shell:
        """
        diamond makedb --in {input.fasta} --db {output.index} --quiet 2> {log}
        """


rule build_dna_taxonomy:
    input:
        vmr_file = vmr_file,
    output:
        dna_tax = os.path.join(out_dir, '{basename}.taxonomy.csv')
    run:
        vmr = pd.read_csv(input.vmr_file, sep = '\t')
        vmr = vmr.rename(columns={'GenBank Assembly ID':'ident', 'Virus name(s)': 'name', 'Exemplar or additional isolate': 'exemplar_or_additional'})
        print(vmr.shape)
        vmr = vmr.dropna(subset=['ident'])
        vmr = vmr.drop_duplicates(subset=['ident']) #only keep first
        print(vmr.shape)
        vmr.set_index('ident', inplace=True)
        vmr['superkingdom'] = 'Viruses'
        tax_columns = ['superkingdom', 'Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', \
                    'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', \
                    'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species', \
                    'exemplar_or_additional', 'name']
        tax_info = vmr[tax_columns]
        tax_info = tax_info.rename(columns=str.lower) # lowercase all names
        tax_info.to_csv(output.dna_tax)


rule build_prot_taxonomy:
    input:
        fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
        fastas=ancient(Checkpoint_MakePattern("{fn}")),
        taxonomy = os.path.join(out_dir, '{basename}.taxonomy.csv'),
    output:
        prot_tax = os.path.join(out_dir, '{basename}.protein-taxonomy.csv'), #columns prot_name, dna_acc, full_lineage
    run:
        # open fromfile
        import screed
        
        # build dna acc --> lineage dict
        dna_acc2lineage = {}
        with open(str(input.taxonomy)) as csvfile:
            r = csv.DictReader(csvfile)  # fields  ['ident', 'superkingdom', 'realm', 'subrealm', 'kingdom', 'subkingdom', 'phylum', \
                                                #    'subphylum', 'class', 'subclass', 'order', 'suborder', \
                                                #    'family', 'subfamily', 'genus', 'subgenus', 'species', \
                                                #    'exemplar_or_additional', 'name']
            for row in r:
                # rename ident to dna_acc
                row['dna_acc'] = row['ident']
                row.pop('ident')
                dna_acc2lineage[row['dna_acc']] = row

        with open(str(output.prot_tax), "a") as outF:
            header = ['ident', 'protein_name', 'dna_acc', 'superkingdom', \
                      'realm', 'subrealm', 'kingdom', 'subkingdom', 'phylum', \
                      'subphylum', 'class', 'subclass', 'order', 'suborder', \
                      'family', 'subfamily', 'genus', 'subgenus', 'species', \
                      'exemplar_or_additional', 'name']
            # open dictwriter and write header
            w = csv.DictWriter(outF, fieldnames=header)
            # grab info from each protein file
            with open(str(input.fromfile)) as csvfile:
                r = csv.DictReader(csvfile) # fields: ["name", "genome_filename", "protein_filename"]
                for row in r:
                    # get the dna assembly accession
                    if row["protein_filename"]:
                        dna_acc = row['name'].split(' ')[0]
                        lineageInfo = dna_acc2lineage[dna_acc]
                        # open protein file with screed and grab the name of each fasta entry
                        with screed.open(str(row["protein_filename"])) as protF:
                            for record in protF:
                                prot_name = record.name
                                # get the prot accession from the protein name
                                prot_acc = prot_name.split(" ")[0]
                                # write to output file
                                w.writerow({"ident": prot_acc, "protein_name": prot_name, **lineageInfo})

