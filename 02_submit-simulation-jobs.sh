# sbatch --array=1-1000 -J neutral 02_slurm-simulate.sh --configfile ./02_simulation-parameters/neutral.yaml

for i in {1..5}
do
    sbatch --array=1-1000 -J hard${i} 02_slurm-simulate.sh --configfile ./02_simulation-parameters/hard.yaml
done

for i in {1..15}
do
    sbatch --array=1-1000 -J rnm${i} 02_slurm-simulate.sh --configfile ./02_simulation-parameters/rnm.yaml
done

for i in {1..15}
do
    sbatch --array=1-1000 -J sgv${i} 02_slurm-simulate.sh --configfile ./02_simulation-parameters/sgv.yaml
done
