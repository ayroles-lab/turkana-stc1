
configfile: '01_config.yaml'

rule all:
    input: "output/sfs/sfs.tsv"

rule concatenate_sfs:
    input:
        expand(
            "raw-data/20220112_raw-sfs/chr{chrom}.sfs.txt",
            chrom=[str(i) for i in range(1, 23)]
        )
    output: "output/sfs/sfs.tsv"
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/prepare-data/concatenate-sfs.tsv"
    log: "notebook-logs/prepare-data/concatenate-sfs.py.ipynb"
    notebook: "notebooks/prepare-data/contatenate-sfs.py.ipynb"
