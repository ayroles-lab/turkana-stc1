
configfile: "03_config.yaml"

rule all:
    input:
        "fig/sfs.pdf",
        "fig/recombination-rates.pdf"

rule plot_recombination_rates:
    input:
        sweep = "output/empirical-statistics/recombination-at-sweep.tsv",
        chrom = "output/empirical-statistics/recombination-at-chromosome-8.tsv"
    output:
        "fig/recombination-rates.pdf"
    conda: "envs/r.yaml"
    notebook: "notebooks/plot/recombination-rates.r.ipynb"

rule plot_sfs:
    input:
        empirical = "output/empirical-statistics/sfs.tsv",
        simulated = "output/simulation-data/simulated-neutral-sfs.tsv"
    output: "fig/sfs.pdf"
    conda: "envs/r.yaml"
    notebook: "notebooks/plot/plot-sfs.r.ipynb"

rule prepare_simulated_sfs:
    output: "output/simulation-data/simulated-neutral-sfs.tsv"
    conda: "envs/simulate.yaml"
    params:
        num_neutral_simulations_to_use = 20
    notebook: "notebooks/plot/prepare-simulated-sfs.py.ipynb"
