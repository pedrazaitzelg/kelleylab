#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=go_analysis                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=aanakamo@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=7                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=15G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

### for running GO analysis in R for each species/tissue

subset=${1}

cd /hb/scratch/aanakamo/kelleylab_rotation/go_analysis

module load R

Rscript ~/kelleylab_rotation/scripts/go_analysis.R -s ${subset}

