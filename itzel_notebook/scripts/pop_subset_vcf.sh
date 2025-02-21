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

## analyses by population

# set working directory
cd /hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/subset_vcf

# set variables
#snpeff vcf file
vcf_file=/hb/groups/kelley_training/itzel/population_bears_proj24/snpEff_out/new_vcf/new_all_genes_ann.vcf
#directory of individual_txt files
dir=/hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/indiv_pop
#specific location file containing individual ids
location=(echo locations.txt | awk ' {print $1; } ')
loc_fil=${dir}/${location}_indiv.txt


#load module
module load vcftools

# include only individuals in file
vcftools ${vcf_file} --keep ${loc_fil} --out ${location}_subset

#unload module
module unload vcftools


