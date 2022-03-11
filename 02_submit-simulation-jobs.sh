sbatch --array=1-1000 -J basic-neutral 02_slurm-simulate.sh --configfile ./02_simulation-parameters/neutral.yaml
sbatch --array=1-1000 -J codominant-hard 02_slurm-simulate.sh --configfile ./02_simulation-parameters/hard.yaml
sbatch --array=1-1000 -J codominant-rnm 02_slurm-simulate.sh --configfile ./02_simulation-parameters/rnm.yaml
sbatch --array=1-1000 -J codominant-sgv 02_slurm-simulate.sh --configfile ./02_simulation-parameters/sgv.yaml
sbatch --array=1-1000 -J dominant-hard 02_slurm-simulate.sh --configfile ./02_simulation-parameters/dominant-hard.yaml
sbatch --array=1-1000 -J dominant-rnm 02_slurm-simulate.sh --configfile ./02_simulation-parameters/dominant-rnm.yaml
sbatch --array=1-1000 -J dominant-sgv 02_slurm-simulate.sh --configfile ./02_simulation-parameters/dominant-sgv.yaml
