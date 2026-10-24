#!/bin/bash
#module load cuda/11.2.2
#module load cudnn/8.1.0.77

test=$1

if [ "$test" = "dry" ]; then
    snakemake -n  --configfile inputs.yaml --cores 32 --jobs 1 --latency-wait 10
else
    snakemake  --configfile inputs.yaml --cores 32 --jobs 1 --latency-wait 10
fi

