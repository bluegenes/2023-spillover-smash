import pandas as pd
import csv

out_dir = "output.spillover"
logs_dir = os.path.join(out_dir, 'logs')
basename="spillover"
db_basename = "vmr_MSL38_v1"

sp_fromfile = 'output.spillover/spillover.fromfile.csv'
sp = pd.read_csv(sp_fromfile)

params = {"dna": {"ksize": [9,15,21], "scaled": 1, "threshold_bp": 0},
          "protein": {"ksize": [5,7], "scaled": 1, "threshold_bp": 0}} # 6,10

TAX_FILE = 'inputs/vmr_MSL38_v1.taxonomy.csv'

onstart:
    print("------------------------------")
    print("sourmash taxonomic classification workflow")
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")

#build parameter combinations for easy expansion below
dna_params = expand(f"{basename}-x-{db_basename}.dna.k{{k}}-sc{{scaled}}.t{{thresh}}", k=params['dna']['ksize'], 
                                                                                       scaled=params['dna']['scaled'],
                                                                                       thresh=params['dna']['threshold_bp'])
prot_params = expand(f"{basename}-x-{db_basename}.protein.k{{k}}-sc{{scaled}}.t{{thresh}}", k=params['protein']['ksize'], 
                                                                                            scaled=params['protein']['scaled'],
                                                                                            thresh=params['protein']['threshold_bp'])
all_params = prot_params
# all_params = dna_params + prot_params


rule all:
    input:
        expand(f"{out_dir}/fastmultigather/{{search_params}}.gather.csv", search_params = all_params), 
        # expand(f"{out_dir}/fastmultigather/{{search_params}}.gather.with-lineages.csv", search_params = all_params), 


rule index_database:
    input:
        db_zip = "output.vmr/{db_basename}.{moltype}.zip"
    output:
        db_rocksdb = "output.vmr/{db_basename}.{moltype}-k{ksize}.rocksdb"
    conda: "conf/env/branchwater.yml"
    log: f"{logs_dir}/index/{{db_basename}}.{{moltype}}-k{{ksize}}.index.log"
    threads: 1
    shell:
        """
        sourmash scripts index {input.db_zip} -m {wildcards.moltype} \
                               -k {wildcards.ksize} --scaled 1 -o {output.db_rocksdb} 2> {log}
        """

rule sourmash_fastmultigather:
    input:
        query_zip = f"{out_dir}/{{basename}}.{{moltype}}.zip",
        database  = "output.vmr/{db_basename}.{moltype}-k{ksize}.rocksdb",
    output:
        f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        time=60,
        partition="low2",
    threads: 64
    log: f"{logs_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    benchmark: f"{logs_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    conda: "conf/env/branchwater.yml"
    shell:
        """
        echo "DB: {input.database}"
        echo "DB: {input.database}" > {log}

        sourmash scripts fastmultigather --threshold {wildcards.threshold} \
                          --moltype {wildcards.moltype} --ksize {wildcards.ksize} \
                          --scaled {wildcards.scaled} {input.query_zip} \
                          {input.database} -o {output} 2> {log}
        """


rule tax_annotate:
    input:
        gather_csv = f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
        lineages = TAX_FILE,
    output:
        f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.with-lineages.csv",
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
        sourmash tax annotate -g {input.gather_csv} -t {input.lineages} --ictv 2> {log}
        """