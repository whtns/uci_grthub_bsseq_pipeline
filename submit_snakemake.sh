#!/bin/bash
#SBATCH --job-name=snakemake_rnaseq
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --error=snakemake-%j.err
#SBATCH --output=snakemake-%j.out

# Load required modules
module load python/3.8

# Activate your conda environment if you have one with snakemake
# conda activate snakemake_env

# Run Snakemake with SLURM integration
snakemake \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem={cluster.mem} --time={cluster.time} --job-name={cluster.job-name} --output={cluster.output} --error={cluster.error}" \
    --cluster-config cluster.yaml \
    --jobs 10 \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds \
    --reason

# Alternative: Use snakemake with slurm profile if available
# snakemake --profile slurm --jobs 10
