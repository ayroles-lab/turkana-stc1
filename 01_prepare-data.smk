
configfile: '03_config.yaml'


rule all:
    input:
        expand("output/empirical-statistics/stc1-{kb}kb_max-pi.txt",
               kb=[1, 10, 100]),
        "output/empirical-statistics/simulation-parameter-estimates.txt",
        "output/empirical-windows/data.tar",
        "output/empirical-windows/logdata.tar",
        "output/empirical-statistics/recombination-at-sweep.tsv",
        "output/empirical-statistics/recombination-at-chromosome-8.tsv",

rule recombination_rates:
    input: "raw-data/20220216_recombination-maps/maps_chr.8"
    output:
        recombination_at_sweep = "output/empirical-statistics/recombination-at-sweep.tsv",
        chromosome_recombinations = "output/empirical-statistics/recombination-at-chromosome-8.tsv"
    conda: "envs/simulate.yaml"
    notebook: "notebooks/prepare-data/recombination.py.ipynb"

rule compress_empirical_log_features:
    input: "output/empirical-windows/npy-log-scale/sweep.npy"
    output: 'output/empirical-windows/logdata.tar'
    params:
        npy_dir = 'output/empirical-windows/npy-log-scale'
    shell:
        "tar -czf {output} {params.npy_dir} ;"


rule compress_empirical_features:
    input: "output/empirical-windows/npy/sweep.npy"
    output: 'output/empirical-windows/data.tar'
    params:
        npy_dir = 'output/empirical-windows/npy'
    shell:
        "tar -czf {output} {params.npy_dir} ;"


rule empirical_window_features:
    input:
        ms = 'output/empirical-windows/ms/{window}.ms',
        normalization_stats = config['stats_file_location']
    output:
        npy = 'output/empirical-windows/npy/{window}.npy',
        log_npy = 'output/empirical-windows/npy-log-scale/{window}.npy',
        features = 'output/empirical-windows/features/{window}.tsv',
        stats = 'output/empirical-windows/features/{window}-stats.tsv',
    params:
        first_position_in_vcf = 23350029,
        outdir = 'output/empirical-windows/npy'
    conda: 'envs/simulate.yaml'
    notebook: 'notebooks/prepare-data/empirical-window-features.py.ipynb'


rule empirical_window_ms:
    input:
        data = config["raw_sweep_region_vcf"]
    output:
        ms = 'output/empirical-windows/ms/sweep.ms'
    conda: 'envs/simulate.yaml'
    script: 'scripts/prepare-data/vcf-to-ms.py'


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
