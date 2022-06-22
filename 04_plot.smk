
configfile: "03_config.yaml"

def relate_sites_of_interest():
    with open("output/clues/05_sample-branch-lengths/stc1-sites-of-interest.txt") as f:
        return [int(line.strip()) for line in f]


rule all:
    input:
        "fig/sfs.pdf",
        "fig/recombination-rates.pdf",
        "fig/sweepfinder.pdf",
        "fig/selection-scan.pdf",
        "fig/sweep-signature_dominant.pdf",
        "fig/sweep-signature_codominant.pdf",
        "fig/clues-results.pdf",
        "output/estimated-sweep-time.txt",
        "fig/sim-timing.pdf",
        "fig/learning-curves.pdf",
        "fig/selstrength-validation.pdf",
        "fig/sweep-mode_validation.pdf"


rule plot_ml_results:
    input:
        sim_params = "output/simulation-data-processed/parameters/codominant_parameters-clean.tsv",
        selstrength = "output/inferences-training/log-sel-strength_codominant_validation_replicate-0.tsv",
        sweepmode = "output/inferences-training/sweep-mode_codominant_validation_replicate-0.tsv"
    output:
        timings_table = "output/estimated-sweep-time.txt",
        timing = "fig/sim-timing.pdf",
        learning_curves = "fig/learning-curves.pdf",
        selstrength = "fig/selstrength-validation.pdf",
        sweepmode = "fig/sweep-mode_validation.pdf"
    params:
        predicted_s = 0.041,
        predicted_f0 = 0.06
    conda: "envs/r.yaml"
    notebook: "notebooks/plot/ml-results.r.ipynb"

rule plot_clues:
    input:
        clues = "output/clues/clues-results-tidy.tsv",
        arg_info = "output/clues/04_popsize-inference/stc1-popsizes.mut"
    output:
        results = "fig/clues-results.pdf",
        individual_plots = expand(
            "fig/clues-individual-sites/clues_{site}.pdf",
            site=relate_sites_of_interest()
        )
    params:
        plot_dir = directory("fig/clues-individual-sites/")
    conda: "envs/r.yaml"
    notebook: "notebooks/plot/clues.r.ipynb"


rule prepare_clues:
    input:
        expand("output/clues/06_clues-results/clues-result_{site}{suffix}",
               site=relate_sites_of_interest(),
               suffix=[".log", ".epochs.npy", ".freqs.npy", ".post.npy"])
    output: "output/clues/clues-results-tidy.tsv"
    conda: "envs/simulate.yaml"
    notebook: "notebooks/plot/prepare-clues-results.py.ipynb"

rule plot_sweep_signals:
    input: "output/sweep-signals/sweep-signals_{datatype}.tsv"
    output: "fig/sweep-signature_{datatype}.pdf",
    conda: "envs/r.yaml"
    notebook: "notebooks/plot/sweep-signals.r.ipynb"

rule prepare_sweep_signals:
    input:
        empirical = "output/empirical-windows/data.tar",
        simulated = "output/simulation-data/{datatype}/data.tar",
        parameters = "output/simulation-data-processed/parameters/{datatype}_parameters-clean.tsv"
    output: "output/sweep-signals/sweep-signals_{datatype}.tsv"
    params:
        num_sel_bins = 3
    conda: "envs/simulate.yaml"
    notebook: "notebooks/prepare-data/prepare-sweep-signals.py.ipynb"

rule plot_selection_scan:
    input: "output/selection-scan/selection-scan-features.tsv"
    output: "fig/selection-scan.pdf"
    conda: "envs/r.yaml"
    notebook: "notebooks/plot/selection-scan.r.ipynb"

rule plot_sweepfinder:
    input: "output/sweepfinder/sweed-results.tsv"
    output: "fig/sweepfinder.pdf"
    conda: "envs/r.yaml"
    params:
        popsize = 30000
    notebook: "notebooks/plot/sweepfinder-results.r.ipynb"

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
