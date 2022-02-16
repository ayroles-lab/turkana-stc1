sbatch --array=1-1000 -J neutral 02_slurm-simulate.sh --configfile ./02_simulation-parameters/neutral.yaml
sbatch --array=1-10000 -J hard 02_slurm-simulate.sh --configfile ./02_simulation-parameters/hard.yaml
sbatch --array=1-20000 -J rnm 02_slurm-simulate.sh --configfile ./02_simulation-parameters/rnm.yaml
sbatch --array=1-20000 -J sgv 02_slurm-simulate.sh --configfile ./02_simulation-parameters/sgv.yaml
