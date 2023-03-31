import pandas as pd

out_dir = "output.vmr"
logs_dir = os.path.join(out_dir, 'logs')

gb_file = 'inputs/VMR_21-221122_MSL37.genbank.taxonomy.csv'
rf_file = 'inputs/VMR_21-221122_MSL37.refseq.taxonomy.csv'

gb = pd.read_csv(gb_file)
rf = pd.read_csv(rf_file)

gb_acc = gb['accession'].tolist()
rf_acc = rf['accession'].tolist()

ACCESSIONS = {'genbank': gb_acc, 'refseq': rf_acc}
all_accessions = gb_acc + rf_acc

databases= ['genbank', 'refseq']
db_accs = []
for db in databases:
    accs = expand("{db}/{acc}", db=db, acc=ACCESSIONS[db])
    db_accs +=accs


rule all:
    input: 
        genomes=expand("{dbacc}.fna", dbacc=db_accs),
        #proteomes=lambda w: expand("{db}/{acc}.faa", db=databases,acc= ACCESSIONS[w.db]),
        #"{db}.dna-sc1.zip",
        #"{db}.protein-sc1.zip"


rule download_dna:
    output: "{db}/{acc}.fna"
    params:
        url = lambda w: "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id={w.acc}&extrafeat=null&conwithfeat=on&&retmode=html&tool=portal&withmarkup=on&&maxdownloadsize=1000000"
        #url = lambda w: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={w.acc}&rettype=fasta&retmode=text",
    log: os.path.join(logs_dir, "downloads", "{db}/{acc}.dna.log")
    benchmark: os.path.join(logs_dir, "downloads", "{db}/{acc}.dna.benchmark")
    threads: 1
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="bml",#low2
    shell:
        """
        curl -o {output} -L {params.url} 2> {log}
        """

rule download_protein:
    output: "{db}/{acc}.faa"
    params:
        url = lambda w: "curl https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={w.acc}&rettype=fasta&retmode=text"
    log: os.path.join(logs_dir, "downloads", "{db}/{acc}.protein.log")
    benchmark: os.path.join(logs_dir, "downloads", "{db}/{acc}.protein.benchmark")
    threads: 1
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="low2",
    shell:
        """
        curl {params.url} -o {output} 2> {log}
        """


# to do: make function instead, since some protein files may not actually exist??
#rule write_fromfile:
#    input: 
#        genomes=lambda w: expand("{db}/{acc}.fna", acc= ACCESSIONS[w.db]),
#        proteomes=lambda w: expand("{db}/{acc}.faa", acc= ACCESSIONS[w.db]),
#    output: "{db}.fromfile.csv"
#    run:
#        with open(str(output), 'w') as out:
#            header = 'ident,genome_file,protein_file'
#            out.write(','.join(header) + '\n')
#            #for inf in input.genomes:
#

