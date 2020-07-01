module load WiggleTools deeptools UCSCtools

read workdir config params <<< "$@"

snakemake --snakefile Snakefile -d $workdir --configfile $config $params
