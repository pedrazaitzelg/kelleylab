#!/bin/bash
#SBATCH --job-name=build_snp
#SBATCH --partition=256x44
#SBATCH --mail-type=ALL
#SBATCH --mail-user=igpedraz@ucsc.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --output=snpbuild.out
#SBATCH --error=snpbuild.err
#SBATCH --mem=20G

## create reference genome in SnpEFF

#directory
cd /hb/home/igpedraz/.conda/envs/snpeff_env/share/snpeff-5.2-1
#program
module load miniconda3
#activate environment
conda activate snpeff_env
#command
#-Xmx200m used to expand memory
snpEff build -gtf22 -v -d -Xmx2000m UrsArc2.0
