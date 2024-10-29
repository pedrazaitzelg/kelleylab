#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=2-00:00:00                # Max time for job to run
#SBATCH --job-name=sbedtool_extract           # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=15G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-35]                   # array job

### script filter vcf files for only biallele sites ###
### bedtools had an issue - "NonRef" ###

#working directory
cd /hb/groups/kelley_training/itzel/fall24/filter_out

#load module
module load vcftools

#set variables
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /hb/groups/kelley_training/itzel/fall24/sorted.txt)
vcf_name=$(echo ${LINE} | awk '{ print $1; }')
vcf_file=/hb/groups/kelley_training/itzel/fall24/$vcf_name

./vcftools --gzvcf $vcf_file --min-alleles 2 --max-alleles 2 --recode --out filtered_out.vcf
