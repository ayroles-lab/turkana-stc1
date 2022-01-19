
CONFIGFILES=$(ls 02_simulation-parameters/*.yaml)

for file in ${CONFIGFILES}
do
    snakemake -c4 --use-conda --snakefile 02_simulate.smk --configfile ${file} "$@" --config slim=bin/slim3.7 normalization_stats=resources/normalization-stats.tsv simulations=3 use_subdirectory=true
done
