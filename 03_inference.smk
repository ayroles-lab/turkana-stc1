
configfile: '03_config.yaml'


rule all:
    input:
        data_reports = expand(
            "output/simulation-data-processed/info/{sim_id}_data-report.txt",
            sim_id=(config["training_ids"])
        ),
        training_inferences = expand(
            "output/inferences-training/{target}_{training}_{testing}.tsv",
            target=config["inference_targets"],
            training=config["training_ids"],
            testing=["training", "validation"]
        )


rule fit_model:
    input:
        training = "output/simulation-data-processed/balanced/{target}_{training}_training.tsv",
        validation = "output/simulation-data-processed/balanced/{target}_{training}_validation.tsv",
        data = "output/simulation-data/{training}/data.tar",
        logdata  = "output/simulation-data/{training}/logdata.tar"
    output:
        fit_model = "output/trained-models/{target}_{training}.pth",
        model_object = "output/trained-models/{target}_{training}.pkl",
        model_labels = "output/trained-models/{target}_{training}_labels.txt",
        training_inferences = "output/inferences-training/{target}_{training}_training.tsv",
        validation_inferences = "output/inferences-training/{target}_{training}_validation.tsv",
        fit_report = "output/model-fitting/{target}_{training}_fit.tsv"
    params:
        save_model = True,
        save_inferences = True,
        use_log_data = False,
        epochs = config["epochs_for_model_training"]
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/inference/fit-neural-network_{target}_{training}.tsv"
    log: "notebook-logs/inference/fit-neural-network_{target}_{training}.py.ipynb"
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
    log: "notebook-logs/inference/combine-simulations_{sim}.py.ipynb"
    notebook: "notebooks/inference/combine-simulations.py.ipynb"
