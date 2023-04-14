import pandas as pd
import csv

out_dir = "output.sp"
logs_dir = os.path.join(out_dir, 'logs')
basename="spillover"

sp_fromfile = 'output.sp/spillover.fromfile.csv'
sp = pd.read_csv(sp_fromfile)
dna_names = sp["name"][sp["genome_filename"].notnull()].tolist()
dna_acc = [n.split(' ')[0] for n in dna_names]
prot_names = sp["name"][sp["protein_filename"].notnull()].tolist()
prot_acc = [n.split(' ')[0] for n in prot_names]

# search_databases = config['search_databases'] # must be dictionary
# search_databases = f"{basename}.{{moltype}}.zip"), moltype = ['dna', 'protein']),
ksize = config.get("ksize", [31])
if not isinstance(ksize, list):
    ksize=[ksize]


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
        expand(os.path.join(out_dir, 'gather', 'protein', '{acc}.k{ks}.gather.csv'), ks=ksize, acc=prot_acc),
        expand(os.path.join(out_dir, 'prefetch', 'dna', '{acc}.k{ks}.prefetch.csv'), ks=ksize, acc=dna_acc),
        expand(os.path.join(out_dir, 'prefetch', 'protein', '{acc}.k{ks}.prefetch.csv'), ks=ksize, acc=prot_acc),
    

rule sourmash_prefetch:
    input:
        query_zip=os.path.join(out_dir, f"{basename}.{{moltype}}.sig.zip"),
        database= "output.vmr/ictv.{moltype}.zip"
        # databases = lambda w: search_databases[f"k{w.ksize}"],
    output:
        prefetch_csv=os.path.join(out_dir, 'prefetch', '{moltype}', '{acc}.k{ksize}.prefetch.csv'),
        prefetch_txt=os.path.join(out_dir, 'prefetch', '{moltype}','{acc}.k{ksize}.prefetch.txt'),
    params:
        threshold_bp = 0,
        alpha_cmd = lambda w: "--" + w.moltype,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        time=10000,
        #partition="bmh",
    log: os.path.join(logs_dir, "prefetch", '{moltype}', "{acc}.k{ksize}.prefetch.log")
    benchmark: os.path.join(logs_dir, "prefetch", '{moltype}', "{acc}.k{ksize}.prefetch.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB(s): {input.database}"
        echo "DB(s): {input.database}" > {log}

        sourmash sig grep {wildcards.acc} {input.query_zip} {params.alpha_cmd} \
                          --ksize {wildcards.ksize} | sourmash prefetch - {input.database} \
                          {params.alpha_cmd} --ksize {wildcards.ksize} --threshold-bp {params.threshold_bp} \
                          -o {output.prefetch_csv} > {output.prefetch_txt} 2>> {log}
        
        touch {output.prefetch_txt}
        touch {output.prefetch_csv}
        """

rule sourmash_gather:
    input:
        query_zip=os.path.join(out_dir, f"{basename}.{{moltype}}.sig.zip"),
        database = "output.vmr/ictv.{moltype}.zip"
    output:
        gather_csv=os.path.join(out_dir, 'gather', '{moltype}', '{acc}.k{ksize}.gather.csv'),
        gather_txt=os.path.join(out_dir, 'gather', '{moltype}', '{acc}.k{ksize}.gather.txt'),
    params:
        threshold_bp = 0,
        alpha_cmd = lambda w: "--" + w.moltype,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        time=10000,
        #partition="bmh",
    log: os.path.join(logs_dir, "gather", '{moltype}', "{acc}.k{ksize}.gather.log")
    benchmark: os.path.join(logs_dir, "gather", '{moltype}', "{acc}.k{ksize}.gather.benchmark")
    conda: "conf/env/sourmash.yml"
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