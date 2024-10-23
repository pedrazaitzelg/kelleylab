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
#SBATCH --array=[1-1]                   # array job

#working directory
cd /hb/groups/kelley_training/itzel/fall24/bedtools_out

#load module
module laod bedtools

#set variables
file_loc=hb/groups/kelley_training/itzel/fall24
vcf_file=$file_loc/brownbear_allsites.*.g.vcf.gz

#script 
bedtools intersect -a  $vcf_file \
                   -b genes_INSIG.bed  -wa | insig_out.txt
                    
