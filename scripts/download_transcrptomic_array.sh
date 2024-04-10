#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=getSRA                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=aanakamo@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=1                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-187]                  # array job

### for paralellizing each SRA sample download into a job array

cd /hb/groups/kelley_training/itzel/data/transcriptomic

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p species_tissue_sra_state.txt)
species=$(echo ${LINE} | awk '{ print $1; }')
tissue=$(echo ${LINE} | awk '{ print $2; }')
sra_acc=$(echo ${LINE} | awk '{ print $3; }')
state=$(echo ${LINE} | awk '{ print $4; }')

echo "downloading sra sample: ${sra_acc} (${species}, ${tissue}, ${state})"

prefetch --output-directory ${species}/${tissue} ${sra_acc}
fastq-dump --outdir ${species}/${tissue} --skip-technical --readids --read-filter pass \
    --gzip --dumpbase --split-3 --clip ${species}/${tissue}/${sra_acc}/${sra_acc}.*
rm -r ${species}/${tissue}/${sra_acc}
