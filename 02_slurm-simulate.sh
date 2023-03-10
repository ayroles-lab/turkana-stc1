#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=18000
#SBATCH --account=bscb02
#SBATCH --output="/home/ivc2/slurm-outputs/%x-%j-%a.out"
#SBATCH --partition=regular,long7,long30

echo "Workstation is ${HOSTNAME}, partition is ${SLURM_JOB_PARTITION}."

# Create working directory and the destination folder for results.
WORKDIR=/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}
DATAHOME=/fs/cbsuclarkfs1/storage/ivc2/turkana
RESULTSHOME=${DATAHOME}/raw-simulations/${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}

# Create relevant directory structure
mkdir -p ${WORKDIR}
cd ${WORKDIR}

# Mount storage
/programs/bin/labutils/mount_server cbsuclarkfs1 /storage

echo "Copying analysis scripts."
cp -r ~/turkana/* .
echo "Linking SLiM executable."
mkdir bin
ln -s ~/bin/slim3.7 bin/slim3.7

# We need this for conda environments to work in a script.
echo "Activating conda environment."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate python39

echo "Running simulations."

snakemake -c1 --use-conda --snakefile 02_simulate.smk "$@" --config slim=bin/slim3.7 normalization_stats=resources/normalization-stats.tsv --conda-prefix ~/turkana-conda

echo "Deactivating conda."
conda deactivate

echo "Moving results into storage."
mkdir -p ${RESULTSHOME}/data
mkdir -p ${RESULTSHOME}/logdata
mkdir -p ${RESULTSHOME}/features
mkdir -p ${RESULTSHOME}/parameters
mkdir -p ${RESULTSHOME}/logs
mkdir -p ${RESULTSHOME}/genotypes
mv output/simulations/data.tar ${RESULTSHOME}/data/data_${SLURM_ARRAY_TASK_ID}.tar
mv output/simulations/logdata.tar ${RESULTSHOME}/logdata/logdata_${SLURM_ARRAY_TASK_ID}.tar
mv output/simulations/features.tar.gz ${RESULTSHOME}/features/features_${SLURM_ARRAY_TASK_ID}.tar.gz
mv output/simulations/parameters.tsv ${RESULTSHOME}/parameters/parameters_${SLURM_ARRAY_TASK_ID}.tsv
tar -czf logs.tar.gz logs/simulations
mv logs.tar.gz ${RESULTSHOME}/logs/logs_${SLURM_ARRAY_TASK_ID}.tar.gz
mv output/simulations/genotypes.tar.gz ${RESULTSHOME}/genotypes/genotypes_${SLURM_ARRAY_TASK_ID}.tar.gz

echo "Cleaning up working directory..."
rm -r $WORKDIR
