import numpy as np
import pandas as pd


def read_data(path):
    return pd.read_table(path, dtype={"swept_frequencies": str})


def save_data(df, path):
    df.to_csv(path, index=False, sep="\t")


def add_regime_kind_column(df):
    kinds = {"hard": "hard", "rnm (true)": "soft", "sgv (true)": "soft"}
    results = df.assign(regime_kind=[
        kinds[x]
        if x in kinds
        else None
        for x in df.sweep_mode])
    return results


def add_sweep_age_column(df):
    age = df.slim_generations.copy()
    age.loc[df.sweep_mode == "sgv (true)"] = age - df.sgv_selection_generation
    results = df.assign(sweep_age = age)
    return results


def add_sgv_ages(df):
    result = add_sweep_age_column(df)
    result = result.loc[result.sweep_mode == "sgv (true)"]
    result = result.assign(
        sgv_drift_time = result.slim_generations - result.sweep_age,
    )
    return result


def balance_hard_vs_soft(df, seed=None):
    results = add_regime_kind_column(df)
    results = results.groupby("regime_kind").sample(
        n=results.regime_kind.value_counts().min(), random_state=seed
    )
    return results


def balance_rnm_vs_sgv(df, seed=None):
    results = df.loc[df.sweep_mode != "hard"]
    results = results.groupby("sweep_mode").sample(
        n=df.sweep_mode.value_counts().min(), random_state=seed
    )
    return results


def balance_sgv_f0(df, seed=None):
    return df.loc[df.sweep_mode == "sgv (true)"]


def balance_rnm_num_mutations(df, seed=None):
    return df.loc[df.sweep_mode == "rnm (true)"]


def balance_sweep_age(df, seed=None):
    return add_sweep_age_column(df)


def balance_sgv_ages(df, seed=None):
    return add_sgv_ages(df)


balancing_functions = {
    "log-sel-strength": None,
    "sweep-mode": None,
    "hard-vs-soft": balance_hard_vs_soft,
    "rnm-vs-sgv": balance_rnm_vs_sgv,
    "sweep-age": balance_sweep_age,
    "sgv-f0": balance_sgv_f0,
    "rnm-num-mutations": balance_rnm_num_mutations,
    "sgv-drift-time": balance_sgv_ages,
    "sgv-total-time": balance_sgv_ages,
}

# The actual columns in the DataFrames corresponding to the true labels for each inference target.
target_columns = {
    "log-sel-strength": "log_selection_coefficient",
    "sweep-mode": "sweep_mode",
    "hard-vs-soft": "regime_kind",
    "rnm-vs-sgv": "sweep_mode",
    "sweep-age": "sweep_age",
    "sgv-f0": "actual_frequency_at_selection",
    "rnm-num-mutations": "sample_num_adaptive_copies",
    "sgv-drift-time": "sgv_drift_time",
    "sgv-total-time": "slim_generations",
}
