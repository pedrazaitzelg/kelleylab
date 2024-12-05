#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=agat_run            # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=1                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

# Working directory
cd /hb/groups/kelley_training/itzel/population_bears_proj24/agat_out

# Load module and activate environment
module load miniconda3
conda activate agat_env

# Extract cds regions from gff/gtf file
# Replace infile.gff and infile.fasta with GFF/GTF file locations
agat_sp_extract_sequences.pl -g $infile.gff -f $infile.fasta -t cds
