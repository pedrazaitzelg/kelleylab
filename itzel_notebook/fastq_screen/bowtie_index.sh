#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=fastqc                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=bowtie2.out            # Standard output and error log
#SBATCH --error=bowtie2.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1]                  # array job

module load bowtie/bowtie2-2.3.2

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p hb/groups/kelley_training/itzel/data/fastq_screen/species_gcf.txt)
species=$(echo ${LINE} | awk '{ print $1; }')
GCF=$(echo ${LINE} | awk '{ print $2; }')

genome_dir=/hb/groups/kelley_training/itzel/genomic/hibernation/${species}
fna=/hb/groups/kelley_training/itzel/genomic/hibernation/${species}/GCF_*
echo "running bowtie index for: ${species} ${GCF}"


bowtie2-build genomes/${species}/GCF_* ${species}
