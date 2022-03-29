
configfile: '03_config.yaml'

def relate_sites_of_interest():
    with open("output/inferences-s-other-methods/clues/stc1-sites-of-interest.txt") as f:
        return [int(line.strip()) for line in f]

rule all:
    input:
        data_reports = expand(
            "output/simulation-data-processed/info/{sim_id}_data-report.txt",
            sim_id=(config["training_ids"])
        ),
        training_inferences_replicates = expand(
            "output/inferences-training/{target}_{training}_{testing}_replicate-{k}.tsv",
            target=config["inference_targets"],
            training=config["training_ids"],
            testing=["training", "validation"],
            k=range(config["num_model_replicates"])
        ),
        empirical_inference_replicates = expand(
            "output/inferences-empirical/summaries/{target}_{training}.tsv",
            target=config["inference_targets"],
            training=config["training_ids"]
        ),
        s_estimate = "output/inferences-s-other-methods/messerneher2012-estimate.txt",
        s_esimate_notebook = "output/inferences-s-other-methods/messerneher2012.html",
        sweepfinder = "output/inferences-s-other-methods/sweepfinder2-results.tsv",
        selection_scan = "output/selection-scan/selection-scan-features.tsv",
        # TEMP
        brlens = expand(
            "output/inferences-s-other-methods/clues/branch-lengths/relate-brlens_{site}.timeb",
            site=[23886711]
        )


rule relate_sample_branch_lengths:
    input:
        anc = "output/inferences-s-other-methods/clues/stc1-popsizes.anc",
        mut = "output/inferences-s-other-methods/clues/stc1-popsizes.mut",
        coal = "output/inferences-s-other-methods/clues/stc1-popsizes.coal"
    output: "output/inferences-s-other-methods/clues/branch-lengths/relate-brlens_{site}.timeb"
    params:
        in_prefix = "output/inferences-s-other-methods/clues/stc1-popsizes",
        out_prefix = "output/inferences-s-other-methods/clues/branch-lengths/relate-brlens_{site}",
        mut_rate = 1.083e-8,
        num_samples = 3
    log: "output/inferences-s-other-methods/clues/branch-lengths/logs/relate-brlens_{site}.log"
    shell:
        "bin/relate/scripts/SampleBranchLengths/SampleBranchLengths.sh "
        "-i {params.in_prefix} "
        "-o {params.out_prefix} "
        "-m {params.mut_rate} "
        "--coal {input.coal} "
        "--format b "
        "--num_samples {params.num_samples} "
        "--first_bp {wildcards.site} "
        "--last_bp {wildcards.site} "
        "--seed 13 &> {log}"


rule relate_estimate_popsize:
    input:
        anc = "output/inferences-s-other-methods/clues/stc1-relate.anc",
        mut = "output/inferences-s-other-methods/clues/stc1-relate.mut",
        poplabels = "output/inferences-s-other-methods/clues/stc1-prepared.poplabels"
    output:
        anc = "output/inferences-s-other-methods/clues/stc1-popsizes.anc",
        mut = "output/inferences-s-other-methods/clues/stc1-popsizes.mut",
        coal = "output/inferences-s-other-methods/clues/stc1-popsizes.coal"
    params:
        in_prefix = "output/inferences-s-other-methods/clues/stc1-relate",
        out_prefix = "output/inferences-s-other-methods/clues/stc1-popsizes",
        mut_rate = 1.083e-8
    log: "output/inferences-s-other-methods/clues/stc1-popsizes.log"
    shell:
        "bin/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh "
        "-i {params.in_prefix} "
        "-m {params.mut_rate} "
        "--poplabels {input.poplabels} "
        "--seed 13 "
        "-o {params.out_prefix} &> {log} ; "


rule relate_run_all:
    input:
        haps = "output/inferences-s-other-methods/clues/stc1-prepared.haps",
        sample = "output/inferences-s-other-methods/clues/stc1-prepared.sample",
        rec_map = "output/inferences-s-other-methods/clues/recombination.map"
    output:
        anc = "output/inferences-s-other-methods/clues/stc1-relate.anc",
        mut = "output/inferences-s-other-methods/clues/stc1-relate.mut"
    params:
        out_prefix = "stc1-relate",
        mut_rate = 1.083e-8,
        hap_pop_size = 60000
    log: "output/inferences-s-other-methods/clues/relate-run-all.log"
    shell: "cd output/inferences-s-other-methods/clues ; "
           "../../../bin/relate/bin/Relate "
           "--mode All "
           "-m {params.mut_rate} "
           "-N {params.hap_pop_size} "
           "--haps stc1-prepared.haps "
           "--sample stc1-prepared.sample "
           "--map recombination.map "
           "--seed 13 "
           "-o {params.out_prefix} &> relate-run-all.log"

rule selection_scan:
    input:
        ms = 'output/empirical-windows/ms/sweep.ms'
    output:
        features = 'output/selection-scan/selection-scan-features.tsv',
        stats = 'output/selection-scan/selection-scan-features-stats.tsv'
    params:
        window_size = 20_000,
        window_step = 10_000
    conda: "envs/simulate.yaml"
    notebook: "notebooks/inference/selection-scan.py.ipynb"


rule sweepfinder2:
    input:
        data = "output/inferences-s-other-methods/sweepfinder-data/stc1-sweepfinder.tsv",
        sfs = "output/inferences-s-other-methods/sweepfinder-data/turkana-sfs.tsv"
    output: "output/inferences-s-other-methods/sweepfinder2-results.tsv"
    log: "output/inferences-s-other-methods/sweepfinder-data/sweepfinder2.log"
    params:
        grid = config["sweepfinder_grid_number"]
    shell: "bin/SweepFinder2 -l {params.grid} {input.data} {input.sfs} {output} &> {log}"


rule convert_messer_neher_notebook:
    input: "output/inferences-s-other-methods/messerneher2012.ipynb"
    output: "output/inferences-s-other-methods/messerneher2012.html"
    conda: "envs/simulate.yaml"
    shell: "jupyter nbconvert --to html {input}"


rule infer_s_messer_neher_2012:
    input:
        ms = "output/empirical-windows/ms/sweep.ms",
    output:
        estimate = "output/inferences-s-other-methods/messerneher2012-estimate.txt",
        notebook="output/inferences-s-other-methods/messerneher2012.ipynb"
    params:
        mut_rate = 1.083e-8,
        rec_rate = 1.083e-8
    conda: "envs/simulate.yaml"
    log: notebook="output/inferences-s-other-methods/messerneher2012.ipynb"
    notebook: "notebooks/inference/estimate-s-messer-neher-2012.py.ipynb"


rule replicate_empirical_inferences_table:
    input:
        expand(
            "output/inferences-empirical/{{target}}_{{training}}_empirical_replicate-{k}.tsv",
            k=range(config["num_model_replicates"])
        )
    output:
        replicates_results = "output/inferences-empirical/summaries/{target}_{training}.tsv",
        replicates_statistics = "output/inferences-empirical/summaries/{target}_{training}_statistics.tsv"
    conda:
        "envs/simulate.yaml"
    notebook: "notebooks/inference/empirical-inference-replicates.py.ipynb"


rule apply_replicate_model_to_empirical_data:
    input:
        fit_model = "output/trained-models/{target}_{training}_replicate-{k}.pth",
        model_object = "output/trained-models/{target}_{training}_replicate-{k}.pkl",
        model_labels = "output/trained-models/{target}_{training}_labels_replicate-{k}.txt",
        data = "output/empirical-windows/data.tar",
        logdata = "output/empirical-windows/logdata.tar"
    output:
        inferences = "output/inferences-empirical/{target}_{training}_empirical_replicate-{k}.tsv"
    params:
        application_type = "empirical"
    conda: "envs/ml.yaml"
    notebook: "notebooks/inference/apply-model.py.ipynb"


rule fit_model_replicate:
    input:
        training = "output/simulation-data-processed/balanced/{target}_{training}_training.tsv",
        validation = "output/simulation-data-processed/balanced/{target}_{training}_validation.tsv",
        data = "output/simulation-data/{training}/data.tar",
        logdata  = "output/simulation-data/{training}/logdata.tar"
    output:
        fit_model = "output/trained-models/{target}_{training}_replicate-{k}.pth",
        model_object = "output/trained-models/{target}_{training}_replicate-{k}.pkl",
        model_labels = "output/trained-models/{target}_{training}_labels_replicate-{k}.txt",
        training_inferences = "output/inferences-training/{target}_{training}_training_replicate-{k}.tsv",
        validation_inferences = "output/inferences-training/{target}_{training}_validation_replicate-{k}.tsv",
        fit_report = "output/model-fitting/{target}_{training}_fit_replicate-{k}.tsv"
    params:
        save_model = True,
        save_inferences = True,
        use_log_data = False,
        epochs = config["epochs_for_model_training"]
    conda: "envs/ml.yaml"
    notebook: "notebooks/inference/fit-neural-network.py.ipynb"


rule balance_training_data:
    input:
        training = "output/simulation-data-processed/train-valid-split/{training}_training.tsv",
        validation = "output/simulation-data-processed/train-valid-split/{training}_validation.tsv"
    output:
        balanced_training = "output/simulation-data-processed/balanced/{target}_{training}_training.tsv",
        balanced_validation = "output/simulation-data-processed/balanced/{target}_{training}_validation.tsv"
    params:
        random_seed = 13
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/inference/balance-training-data_{target}_{training}.tsv"
    script: "scripts/inference/balance-training-data.py"

    
rule train_validation_split:
    input:
        sim_params = "output/simulation-data-processed/parameters/{training}_parameters-clean.tsv"
    output:
        training = "output/simulation-data-processed/train-valid-split/{training}_training.tsv",
        validation = "output/simulation-data-processed/train-valid-split/{training}_validation.tsv"
    params:
        random_seed = 13
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/inference/train-validation-split_{training}.tsv"
    script: "scripts/inference/train-validation-split.py"


rule clean_datasets:
    input:
        sim_params = "output/simulation-data/{sim_id}/parameters.tsv"
    output:
        cleaned_parameters = "output/simulation-data-processed/parameters/{sim_id}_parameters-clean.tsv",
        sweep_mode_report = "output/simulation-data-processed/info/{sim_id}_sweep-modes.txt",
        successful_report = "output/simulation-data-processed/info/{sim_id}_data-report.txt",
        failed_report = "output/simulation-data-processed/info/{sim_id}_failed-simulations-report.txt"
    params:
        random_seed = 13
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/inference/clean-dataset_{sim_id}.tsv"
    script: "scripts/inference/clean-dataset.py"


rule combine_simulation_tasks:
    output:
        data = "output/simulation-data/{sim}/data.tar",
        parameters = "output/simulation-data/{sim}/parameters.tsv",
        info = "output/simulation-data/{sim}/info.txt",
        features = "output/simulation-data/{sim}/features.tar.gz",
        logdata  = "output/simulation-data/{sim}/logdata.tar"
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/inference/combine-simulations_{sim}.tsv"
    log: "logs/inference/combine-simulations_{sim}.py.ipynb"
    notebook: "notebooks/inference/combine-simulations.py.ipynb"
