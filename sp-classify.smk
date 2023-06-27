import pandas as pd
import csv

out_dir = "output.spillover"
logs_dir = os.path.join(out_dir, 'logs')
basename="spillover"
#basename="spillover100"

sp_fromfile = 'output.spillover/spillover.fromfile.csv'
#sp_fromfile = 'output.spillover/spillover100.fromfile.csv'
sp = pd.read_csv(sp_fromfile)
dna_names = sp["name"][sp["genome_filename"].notnull()].tolist()
#dna_acc = [n.split(' ')[0] for n in dna_names][:5] # first five for testing
dna_acc = [n.split(' ')[0] for n in dna_names]
prot_names = sp["name"][sp["protein_filename"].notnull()].tolist()
prot_acc = [n.split(' ')[0] for n in prot_names]


# test simian foamy virus
sm_file = 'inputs/simian-foamy.spillover.csv'
sm = pd.read_csv(sm_file)
dna_acc = sm['AccessionNumber']


# search_databases = config['search_databases'] # must be dictionary
# search_databases = f"{basename}.{{moltype}}.zip"), moltype = ['dna', 'protein']),
ksize = config.get("ksize", [21, 31])
if not isinstance(ksize, list):
    ksize=[ksize]

TAX_FILE = 'inputs/VMR_21-221122_MSL37.taxonomy.csv'

onstart:
    print("------------------------------")
    print("sourmash taxonomic classification workflow")
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")


rule all:
    input:
        expand(os.path.join(out_dir, 'gather', 'dna', '{acc}.k{ks}.gather.csv'), ks=ksize, acc=dna_acc),
        expand(os.path.join(out_dir, 'gather', 'dna', '{acc}.k{ks}.gather.classifications.csv'), ks=ksize, acc=dna_acc),
        #expand(os.path.join(out_dir, 'gather', 'protein', '{acc}.k{ks}.gather.csv'), ks=ksize, acc=prot_acc),
        expand(os.path.join(out_dir, 'prefetch', 'dna', '{acc}.k{ks}.prefetch.csv'), ks=ksize, acc=dna_acc),
        #expand(os.path.join(out_dir, 'prefetch', 'protein', '{acc}.k{ks}.prefetch.csv'), ks=ksize, acc=prot_acc),
        expand(os.path.join(out_dir, 'gather', '{basename}.{moltype}-k{ks}.gather-classifications.csv'), basename = basename, moltype='dna', ks=ksize, acc=dna_acc)
    
rule sourmash_prefetch:
    input:
        query_zip=os.path.join(out_dir, f"{basename}.{{moltype}}.zip"),
        database= "output.vmr/ictv.{moltype}.zip"
        # databases = lambda w: search_databases[f"k{w.ksize}"],
    output:
        prefetch_csv=os.path.join(out_dir, 'prefetch', '{moltype}', '{acc}.k{ksize}.prefetch.csv'),
        prefetch_txt=os.path.join(out_dir, 'prefetch', '{moltype}','{acc}.k{ksize}.prefetch.txt'),
    params:
        threshold_bp = 0,
        alpha_cmd = lambda w: "--" + w.moltype,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        time=60,
        partition="low2",
    log: os.path.join(logs_dir, "prefetch", '{moltype}', "{acc}.k{ksize}.prefetch.log")
    benchmark: os.path.join(logs_dir, "prefetch", '{moltype}', "{acc}.k{ksize}.prefetch.benchmark")
    conda: "conf/env/sourmash-dev.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB(s): {input.database}"
        echo "DB(s): {input.database}" > {log}

        sourmash sig grep {wildcards.acc} {input.query_zip} \
                          --ksize {wildcards.ksize} | sourmash prefetch - {input.database} \
                          {params.alpha_cmd} --ksize {wildcards.ksize} --threshold-bp {params.threshold_bp} \
                          -o {output.prefetch_csv} > {output.prefetch_txt} 2>> {log}
        
        touch {output.prefetch_txt}
        touch {output.prefetch_csv}
        """

rule sourmash_gather:
    input:
        query_zip=os.path.join(out_dir, f"{basename}.{{moltype}}.zip"),
        database = "output.vmr/ictv.{moltype}.zip"
    output:
        gather_csv=os.path.join(out_dir, 'gather', '{moltype}', '{acc}.k{ksize}.gather.csv'),
        gather_txt=os.path.join(out_dir, 'gather', '{moltype}', '{acc}.k{ksize}.gather.txt'),
    params:
        threshold_bp = 0,
        alpha_cmd = lambda w: "--" + w.moltype,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        time=240,
        partition="bml",
    log: os.path.join(logs_dir, "gather", '{moltype}', "{acc}.k{ksize}.gather.log")
    benchmark: os.path.join(logs_dir, "gather", '{moltype}', "{acc}.k{ksize}.gather.benchmark")
    conda: "conf/env/sourmash-dev.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB(s): {input.database}"
        echo "DB(s): {input.database}" > {log}

        sourmash sig grep {wildcards.acc} {input.query_zip} {params.alpha_cmd} \
                          --ksize {wildcards.ksize} | sourmash gather - {input.database} \
                          {params.alpha_cmd} --ksize {wildcards.ksize} --threshold-bp {params.threshold_bp} \
                 -o {output.gather_csv} > {output.gather_txt} 2>> {log}
        
        touch {output.gather_txt}
        touch {output.gather_csv}
        """

rule tax_genome:
    input:
        gather_csv=os.path.join(out_dir, 'gather', '{moltype}', '{acc}.k{ksize}.gather.csv'),
        lineages = TAX_FILE,
    output:
        classif_csv=os.path.join(out_dir, 'gather', '{moltype}', '{acc}.k{ksize}.gather.classifications.csv'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, f'gather/{w.moltype}'),
        out_base= lambda w: f'{w.acc}.k{w.ksize}.gather',
    conda: "conf/env/sourmash-ictv.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax genome -g {input.gather_csv} -t {input.lineages} -o {params.out_base} \
                                --output-dir {params.outd} --output-format human csv_summary \
                                --ictv --rank species --force 2> {log}
        touch {output.classif_csv}
        """


rule aggregate_classifications:
    input:
        cl = lambda w: expand(os.path.join(out_dir, 'gather', '{moltype}', '{acc}.k{ksize}.gather.classifications.csv'), acc = dna_acc, moltype = 'dna', ksize=w.ksize),
    output:
        cl_all = os.path.join(out_dir, 'gather', '{basename}.{moltype}-k{ksize}.gather-classifications.csv'),
    run:
        with open(str(output), 'w') as outF:
            w = csv.writer(outF)
            header = None
            empties = []
            for inF in input:
                inF = str(inF)
                if os.stat(inF).st_size == 0: # file is empty
                    acc = inF.split('.', 1)[0]
                    signame = sp[sp["acc"] == acc]["name"]
                    status = "nomatch"
                    empties.append([signame, status])
                with open(inF, "r", newline="") as cl:
                    r = csv.reader(cl)
                    if header is None:
                        header = next(r)
                        w.writerow(header)
                    else:
                        next(r)  # Skip header row
                    for row in r:
                        w.writerow(row)
            for e in empties:
                w.writerow(e)
