NUM_JOBS=1000

for file in $(ls ../simulation-parameters/*.yaml)
do
    NAME=$(basename ${file} .yaml)
    CMD= "sbatch --array=1-${NUM_JOBS} -J ${NAME} slurm-simulate.sh --configfile ${file}"
    echo ${CMD}
    # sbatch --array=1-${NUM_TRAINING_JOBS} -J ${NAME} slurm-simulate.sh --configfile ${file}
done
