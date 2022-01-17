
configfile: '01_config.yaml'


rule all:
    input: "fig/data/sfs.pdf"


##################
# PLOTTING RULES #
##################

rule plot_sfs:
    input: "output/sfs/sfs.tsv"
    output: "fig/data/sfs.pdf"
    conda: "envs/r.yaml"
    notebook: "notebooks/plot/plot-sfs.r.ipynb"
    

##########################
# DATA PREPARATION RULES # 
##########################

rule concatenate_sfs:
    input:
        expand(
            "raw-data/20220112_raw-sfs/chr{chrom}.sfs.txt",
            chrom=[str(i) for i in range(1, 23)]
        )
    output: "output/sfs/sfs.tsv"
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/prepare-data/concatenate-sfs.tsv"
    notebook: "notebooks/prepare-data/contatenate-sfs.py.ipynb"
