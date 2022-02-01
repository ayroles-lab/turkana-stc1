
configfile: '01_config.yaml'


rule all:
    input:
        expand("output/empirical-statistics/stc1-{kb}kb_max-pi.txt",
               kb=[1, 10, 100]),
        "output/empirical-statistics/simulation-parameter-estimates.txt"


rule estimate_simulation_parameters:
    input: "output/empirical-statistics/sfs.tsv"
    output: "output/empirical-statistics/simulation-parameter-estimates.txt"
    params:
        total_genome_size = 3e9,
        reference_ne = 30_000,
        mu_over_r = 1,
        sample_size = 220
    conda: "envs/simulate.yaml"
    notebook: "notebooks/prepare-data/estimate-simulation-parameters.py.ipynb"


rule empirical_statistics:
    input: "raw-data/stc1.vcf.gz"
    output:
        "output/empirical-statistics/stc1-{winsize}kb_max-pi.txt"
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/prepare-data/windowed-pi-{winsize}kb.tsv"
    shell:
        'vcftools --gzvcf {input} --window-pi {wildcards.winsize}000 --window-pi-step 1000 --out "output/empirical-statistics/stc1-{wildcards.winsize}kb" ;'
        "cd output/empirical-statistics ;"
        "sed '1d' stc1-{wildcards.winsize}kb.windowed.pi | sort -gr -k 5 > stc1-{wildcards.winsize}kb_max-pi.txt"

        
rule compress_raw_data:
    input: config["raw_sweep_region_vcf"]
    output: "raw-data/stc1.vcf.gz"
    shell:
        "gzip -c {input} > {output};"


rule concatenate_sfs:
    input:
        expand(
            "raw-data/20220112_raw-sfs/chr{chrom}.sfs.txt",
            chrom=[str(i) for i in range(1, 23)]
        )
    output: "output/empirical-statistics/sfs.tsv"
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/prepare-data/concatenate-sfs.tsv"
    notebook: "notebooks/prepare-data/contatenate-sfs.py.ipynb"
