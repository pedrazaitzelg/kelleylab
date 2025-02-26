#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=filter_pop                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL


## Subset provided populations file into each population
## to be used with the larger vcf file for pop-level stats


#Set Variables
bear_in=/hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/mod_bear_pop.txt

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p $bear_in)
indiv_id=$(echo ${LINE} | awk '{ print $1; }')
#variable not needed
#location=$(echo ${LINE} | awk '{ print $2; }')
sub_loc=$(echo ${LINE} | awk '{ print $3; }')

#create temp file containing location names
awk '{print $2}' ${bear_in} | sort | uniq > locations.txt
location_name=$(echo locations.txt | awk '{ print $1; }')

# create subset files
while read location; do
    # Create a subset for each unique location
    awk -v location="$location" '$2 == location' ${bear_in} | awk '{print $1}' > "${location}_indiv.txt"
done < locations.txt

#remove temp file
rm locations.txt
