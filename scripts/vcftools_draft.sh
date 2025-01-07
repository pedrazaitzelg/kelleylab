#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=vcftools_allele                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

#set variables
ANNOT_VCF=/hb/groups/kelley_training/itzel/population_bears_proj24/vcftool_out/all_genes_ann.vcf
OUT=/hb/groups/kelley_training/itzel/population_bears_proj24/vcftool_out

#load vcftools
module load vcftools

#command for calculating allele frequency for each variant from annotated vcf file
vcftools --gzvcf $ANNOT_VCF --freq2 --out $OUT --max-alleles 2
