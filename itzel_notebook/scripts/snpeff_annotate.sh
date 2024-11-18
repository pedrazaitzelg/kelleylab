#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=2-00:00:00                # Max time for job to run
#SBATCH --job-name=snpeff_run           # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=15G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

### load module ###
module load miniconda3
conda activate snpeff

### locate directory with snpeff location ###
cd /hb/groups/kelley_training/itzel/population_bears_proj24/snpEff_out

### create annotated file destination ###
touch all_genes_ann.vcf

###cd /hb/home/igpedraz/.conda/envs/snpeff

### run vcf file through snpeff to annotate variations ###
snpEff UrsArc2.0 \
/hb/groups/kelley_training/itzel/population_bears_proj24/bedtools_out/allgenes_intersect.vcf \
> /hb/groups/kelley_training/itzel/population_bears_proj24/snpEff_out/all_genes_ann.vcf

### deactivate snpeff ###
conda deactivate
