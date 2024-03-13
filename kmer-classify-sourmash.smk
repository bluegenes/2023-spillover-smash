import pandas as pd
import csv

out_dir = "output.spillover"
logs_dir = os.path.join(out_dir, 'logs')
basename="spillover"
# db_basename = "vmr_MSL38_v1"
# db_dir = "output.vmr"
db_basename = "genbank.2023-05.viral"
db_dir = "output.genbank-viral"

sp_fromfile = 'output.spillover/spillover.fromfile.csv'
sp = pd.read_csv(sp_fromfile)

params = {"dna": {"ksize": [31], "scaled": 1, "threshold_bp": 0}, # 9, 15, 21
          "protein": {"ksize": [5, 7], "scaled": 1, "threshold_bp": 0}} #, 6,10

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
all_params = dna_params # prot_params
# all_params = dna_params + prot_params


rule all:
    input:
        expand(f"{out_dir}/fastmultigather/{{search_params}}.gather.csv", search_params = all_params), 
        expand(f"{out_dir}/fastmultigather/{{search_params}}.gather.with-lineages.csv", search_params = all_params), 


rule index_database:
    input:
        db_zip = "{db_dir}/{db_basename}.{moltype}.zip"
    output:
        db_rocksdb = directory(f"{db_dir}/{{db_basename}}.{{moltype}}-k{{ksize}}.rocksdb"),
        current_file = f"{db_dir}/{{db_basename}}.{{moltype}}-k{{ksize}}.rocksdb/CURRENT",
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
        database  = f"{db_dir}/{{db_basename}}.{{moltype}}-k{{ksize}}.rocksdb/CURRENT",
    output:
        f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        time=60,
        partition="low2",
    threads: 100
    log: f"{logs_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    benchmark: f"{logs_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    conda: "conf/env/branchwater.yml"
    params:
        db_dir = lambda w: f'{db_dir}/{w.db_basename}.{w.moltype}-k{w.ksize}.rocksdb',
    shell:
        """
        echo "DB: {params.db_dir}"
        echo "DB: {params.db_dir}" > {log}

        sourmash scripts fastmultigather --threshold {wildcards.threshold} \
                          --moltype {wildcards.moltype} --ksize {wildcards.ksize} \
                          --scaled {wildcards.scaled} {input.query_zip} \
                          {params.db_dir} -o {output} 2>> {log}
        """


rule tax_annotate_fastmultigather:
    input:
        gather_csv = f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
        lineages = TAX_FILE,
    output:
        f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.with-lineages.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    conda: "conf/env/sourmash-ictv.yml"
    log: f"{logs_dir}/annotate/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.ictv.log"
    params:
        output_dir = f"{out_dir}/fastmultigather",
    shell:
        """
        sourmash tax annotate -g {input.gather_csv} \
                              -t {input.lineages} \
                              --ictv -o {params.output_dir} 2> {log}
        """


rule sourmash_gather:
    input:
        query_zip = ancient(f"{out_dir}/{basename}.{{moltype}}.zip"),
        db_zip = "{db_dir}/{db_basename}.{moltype}.zip",
        fastmultigather =  f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv" ,
    output:
        csv=f"{out_dir}/gather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
        txt=f"{out_dir}/gather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        time=60,
        partition="low2",
    threads: 1
    log: f"{logs_dir}/gather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    benchmark: f"{logs_dir}/gather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    conda: "cbranch",
    params:
        alpha_cmd = lambda w: f"--{w.moltype}",
    shell:
        '''
        echo "DB: {input.db_zip}"
        echo "DB: {input.db_zip}" > {log}
        sourmash sig grep {params.name_regex:q} {input.query_zip} --ksize {wildcards.ksize} {params.alpha_cmd} | \
                          sourmash gather - {input.picklist_zip} \
                          --threshold {wildcards.threshold} \
                          {params.alpha_cmd} --ksize {wildcards.ksize} \
                          --scaled {wildcards.scaled}  \
                           -o {output.csv} > {output.txt} 2>> {log}
        '''
                          #--picklist {input.sample_picklist}:match_md5:md5 \
        # echo "DB: {input.database}"
        # echo "DB: {input.database}" > {log}
        #    echo "DB: {input.picklist_zip}"
        # echo "DB: {input.picklist_zip}" > {log}
        # sourmash sig extract --name "{wildcards.sample}" {input.query_zip} --ksize {wildcards.ksize} {params.alpha_cmd} | \
        #                   sourmash gather - {input.picklist_zip} \
        #                   --threshold {wildcards.threshold} \
        #                   {params.alpha_cmd} --ksize {wildcards.ksize} \
        #                   --scaled {wildcards.scaled}  \
        #                    -o {output} 2>> {log}


## to do: make lineages file for refseq db (or will full genbank one work?)

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
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    log: f"{logs_dir}/annotate/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    benchmark: f"{logs_dir}/annotate/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    params:
        output_dir = f"{out_dir}/fastmultigather",
    shell:
        """
        sourmash tax annotate -g {input.gather_csv} \
                              -t {input.lineages} \
                              -o {params.output_dir} 2> {log}
        """