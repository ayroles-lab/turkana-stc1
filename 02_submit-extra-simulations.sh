for i in {1..5}
do
    sbatch --array=1-1000 -J domhard${i} 02_slurm-simulate.sh --configfile ./02_simulation-parameters/dominant/hard.yaml
done

for i in {1..5}
do
    sbatch --array=1-1000 -J domrnm${i} 02_slurm-simulate.sh --configfile ./02_simulation-parameters/dominant/rnm.yaml
done

for i in {1..5}
do
    sbatch --array=1-1000 -J domsgv${i} 02_slurm-simulate.sh --configfile ./02_simulation-parameters/dominant/sgv.yaml
done

for i in {1..5}
do
    sbatch --array=1-1000 -J widehard${i} 02_slurm-simulate.sh --configfile ./02_simulation-parameters/wide/hard.yaml
done

for i in {1..5}
do
    sbatch --array=1-1000 -J widernm${i} 02_slurm-simulate.sh --configfile ./02_simulation-parameters/wide/rnm.yaml
done

for i in {1..5}
do
    sbatch --array=1-1000 -J widesgv${i} 02_slurm-simulate.sh --configfile ./02_simulation-parameters/wide/sgv.yaml
done
