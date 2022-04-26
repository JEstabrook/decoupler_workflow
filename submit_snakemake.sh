#!/usr/bin/bash
#SBATCH --time 36:00:00
#SBATCH --partition exacloud 
#SBATCH --job-name workflow_submission
#SBATCH --output=logs/workflow_submission_%j.log
#SBATCH --mem 9000

export OMP_NUM_THREADS=1

snakemake -j 125 --use-singularity --rerun-incomplete --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -o {cluster.o} -e {cluster.e} -J {cluster.J} -c {cluster.c} --mem {cluster.mem} --exclude='exanode-4-15,exanode-4-30,exanode-6-27,exanode-6-47,exanode-2-44,exanode-8-[22-24],exanode-3-[0-6],exanode-1-[32-44],exanode-0-[26-29]'" -s Snakefile --latency-wait 1000 --keep-going 
