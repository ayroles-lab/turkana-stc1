NUM_JOBS=500

for file in $(ls ./02_simulation-parameters/*.yaml)
do
    NAME=$(basename ${file} .yaml)
    # CMD="sbatch --array=1-${NUM_JOBS} -J ${NAME} slurm-simulate.sh --configfile ${file}"
    # echo ${CMD}
    sbatch --array=1-${NUM_JOBS} -J ${NAME} 02_slurm-simulate.sh --configfile ${file}
done
