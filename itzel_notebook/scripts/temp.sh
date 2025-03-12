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


## Filter .frq files generated from vcf_tools to remove
## things fixed in both populations an doutput only the differences

#set working directory
cd /hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/allele_freq/temp

##pairwise comparisons

## start with ABC and Alaska

#set variables
#LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p locations.txt)
dir=/hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/indv_stats

## sets location name from locations.txt
location_1=$( awk ' NR==1{print $1} ' locations.txt)
location_2=$( awk ' NR==2{print $1} ' locations.txt)


allele_file1=${dir}/${location_1}/${location_1}.frq  #location of files
allele_file2=${dir}/${location_2}/${location_2}.frq  #location of files

awk ' {print $1; } ' ${allele_file1} | sort | uniq > chrom1.txt  #generate txt file with list of names

chrom=$( awk ' NR==2{print $1} ' chrom1.txt)  #creates variable set to a chrom name

# search file for chromosome and output new file with only that chromosome position
grep ${chrom} ${allele_file1} > ${chrom}_1.txt
grep ${chrom} ${allele_file2} > ${chrom}_2.txt


## sets variable to values in second column
pos=$(echo ${chrom}_1.txt | awk ' {print $2; } ' )
freq=


