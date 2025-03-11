#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=allele_filter                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#sbatch --array[1-2]                     # to run array

## Filter .frq files generated from vcf_tools to remove
## things fixed in both populations an doutput only the differences

#set working directory
cd /hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/allele_freq

##pairwise comparisons

## start with ABC and Alaska

#set variables
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p locations.txt)
location_1=$(echo ${LINE} | awk ' {print $1; 1} ')
location_2=$(echo ${LINE} | awk ' {print $1; 2} ')


loc_file1=$(/hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/indv_stats/${location_1}/${location_1}.frq)
loc_file1=$(/hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/indv_stats/${location_2}/${location_2}.frq)

# generate list of chromosomes


awk ' {print $1; } ' ${loc_file1} | sort | uniq > chrom1.txt
awk ' {print $1; } ' ${loc_file2} | sort | uniq > chrom2.txt

# exclude everything where pop1 - pop2 is 0

for chrom in chrom1.txt
  do


    
