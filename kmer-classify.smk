import pandas as pd
import csv

out_dir = "output.spillover"
logs_dir = os.path.join(out_dir, 'logs')
basename="spillover"
db_basename = "vmr_MSL38_v1"
db_dir = "output.vmr"
# db_basename = "genbank.2023-05.viral"
# db_dir = "output.genbank-viral"

sp_fromfile = 'output.spillover/spillover.fromfile.csv'
sp = pd.read_csv(sp_fromfile)

params = {"dna": {"ksize": [21,31], "scaled": 1, "threshold_bp": 3},
          "protein": {"ksize": [7, 10], "scaled": 1, "threshold_bp": 5}} #, 6,10

TAX_FILE = 'output.vmr/vmr_MSL38_v1.taxonomy.csv'
ANI_THRESHOLD = 0.95
# AAI_THRESHOLD = 0.8 # for protein


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
# all_params = prot_params
# all_params = dna_params
all_params = dna_params + prot_params


rule all:
    input:
        expand(f"{out_dir}/fastmultigather/{{search_params}}.fmgather.csv", search_params = all_params), 
        expand(f"{out_dir}/fastmultigather/{{search_params}}.fmgather.with-lineages.csv", search_params = all_params), 
        expand(f"{out_dir}/fastmultigather/{{search_params}}.classifications.csv", search_params=all_params),
        expand(f"{out_dir}/spillover.{{moltype}}.pairwise.csv", moltype=["dna", "protein"]),
        expand(f"{out_dir}/spillover.dna.cluster.t{{threshold}}.csv", threshold = ANI_THRESHOLD),
        expand(f"{out_dir}/spillover.protein.cluster.t{{threshold}}.csv", threshold = AAI_THRESHOLD),
        expand(f"{out_dir}/spillover.dna.cluster.t{{threshold}}.with-lineages.csv", threshold = ANI_THRESHOLD),
        expand(f"{out_dir}/spillover.protein.cluster.t{{threshold}}.with-lineages.csv", threshold = AAI_THRESHOLD),


rule index_database:
    input:
        db_zip = f"{db_dir}/{{db_basename}}.{{moltype}}.zip"
    output:
        db_rocksdb = directory(f"{db_dir}/{{db_basename}}.{{moltype}}-k{{ksize}}-sc{{scaled}}.rocksdb"),
        current_file = f"{db_dir}/{{db_basename}}.{{moltype}}-k{{ksize}}-sc{{scaled}}.rocksdb/CURRENT",
    conda: "conf/env/branchwater.yml"
    log: f"{logs_dir}/index/{{db_basename}}.{{moltype}}-k{{ksize}}-sc{{scaled}}.index.log"
    benchmark: f"{logs_dir}/index/{{db_basename}}.{{moltype}}-k{{ksize}}-sc{{scaled}}.index.benchmark"
    threads: 1
    params:
        moltype = lambda w: w.moltype.upper() if w.moltype == "dna" else w.moltype,
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
        time=240,
        partition="low2",
    shell:
        """
        sourmash scripts index {input.db_zip} -m {params.moltype} \
                               -k {wildcards.ksize} --scaled {wildcards.scaled} \
                               -o {output.db_rocksdb} 2> {log}
        """

rule sourmash_fastmultigather:
    input:
        query_zip = f"{out_dir}/{{basename}}.{{moltype}}.zip",
        database  = f"{db_dir}/{{db_basename}}.{{moltype}}-k{{ksize}}-sc{{scaled}}.rocksdb/CURRENT",
    output:
        f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.fmgather.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        time=240,
        partition="low2",
    threads: 10
    log: f"{logs_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    benchmark: f"{logs_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    conda: "conf/env/branchwater.yml"
    params:
        db_dir = lambda w: f"{db_dir}/{w.db_basename}.{w.moltype}-k{w.ksize}-sc{w.scaled}.rocksdb",
        moltype = lambda w: w.moltype.upper() if w.moltype == "dna" else w.moltype,
    shell:
        """
        echo "DB: {params.db_dir}"
        echo "DB: {params.db_dir}" > {log}

        sourmash scripts fastmultigather --threshold {wildcards.threshold} \
                          --moltype {params.moltype} --ksize {wildcards.ksize} \
                          --scaled {wildcards.scaled} {input.query_zip} \
                          {params.db_dir} -o {output} 2>> {log}
        """


rule tax_annotate:
    input:
        gather_csv = f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.fmgather.csv",
        lineages = TAX_FILE,
    output:
        f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.fmgather.with-lineages.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    conda: "conf/env/branchwater.yml"
    log: f"{logs_dir}/annotate/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.ictv.log"
    benchmark: f"{logs_dir}/annotate/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.ictv.benchmark"
    params:
        output_dir = f"{out_dir}/fastmultigather",
    shell:
        """
        sourmash tax annotate -g {input.gather_csv} \
                              -t {input.lineages} \
                              --ictv -o {params.output_dir} 2> {log}
        """

rule tax_genome:
    input:
        gather_csv = f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.fmgather.csv",
        lineages = TAX_FILE,
    output:
        f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.classifications.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    conda: "conf/env/sourmash-ictv.yml"
    log: f"{logs_dir}/tax/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.ictv.log"
    benchmark: f"{logs_dir}/tax-genome/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.ictv.benchmark"
    params:
        output_dir = f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}"
    shell:
        """
        sourmash tax genome -g {input.gather_csv} \
                            -t {input.lineages} \
                            --ictv -o {params.output_dir} \
                            --ani-threshold 0.2 --force  2> {log}
        """

rule spillover_pairwise:
    input:
        spillover_zip = f"{out_dir}/{{basename}}.{{moltype}}.zip",
    output:
        f"{out_dir}/spillover.{{moltype}}.pairwise.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        partition="low2",
        time=240,
    conda: "conf/env/branchwater.yml"
    log: f"{logs_dir}/pairwise/spillover.{{moltype}}.pairwise.log"
    benchmark: f"{logs_dir}/pairwise/spillover.{{moltype}}.pairwise.benchmark"
    shell:
        """
        sourmash scripts pairwise {input.spillover_zip} --ani \
                                  --write-all -o {output} 2> {log}
        """

rule spillover_cluster:
    input:
        pairwise: f"{out_dir}/spillover.{{moltype}}.pairwise.csv",
    output:
        f"{out_dir}/spillover.{{moltype}}.cluster.t{threshold}.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        partition="low2",
        time=240,
    conda: "conf/env/branchwater.yml"
    log: f"{logs_dir}/cluster/spillover.{{moltype}}.cluster.t{threshold}.log"
    benchmark: f"{logs_dir}/cluster/spillover.{{moltype}}.cluster.t{threshold}.benchmark"
    shell:
        """
        sourmash scripts cluster {input.pairwise} --threshold {wildcards.threshold} \
                                 -o {output} 2> {log}
        """

rule annotate_clusters:
    input:
        cluster_csv = f"{out_dir}/spillover.{{moltype}}.cluster.t{threshold}.csv",
        lineages = TAX_FILE,
    output:
        f"{out_dir}/spillover.{{moltype}}.cluster.t{threshold}.with-lineages.csv"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    conda: "conf/env/branchwater.yml"
    log: f"{logs_dir}/cluster/spillover.cluster.t{threshold}.annotate.log"
    benchmark: f"{logs_dir}/cluster/spillover.cluster.t{threshold}.annotate.benchmark"
    params:
        output_dir = f"{out_dir}",
    shell:
        """
        python annotate-clusters.py {input.cluster_csv} \
                                    -t {input.lineages} \
                                    --ictv -o {params.output_dir} 2> {log}
        """
