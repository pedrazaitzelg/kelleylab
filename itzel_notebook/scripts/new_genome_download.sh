#!/bin/bash
#SBATCH --job-name=genomic_data
#SBATCH --partition=128x24
#SBATCH --mail-type=ALL
#SBATCH --mail-user=igpedraz@ucsc.edu
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output=genomic_data.out
#SBATCH --error=genomic_data.err
#SBATCH --mem=250M

## Run this script to download the genomic data for each species

## Load module
module load miniconda3

## Activate ncbi_datasets env
conda activate ncbi_datasets

## Options for datasets command
## --include genome,protein,gff3,rna,cds,gtf
## replace $ with corresponding names
datasets download genome accession "${gcf}" --include genome,protein,gff3,rna,cds --filename ${sp}.zip
unzip ${sp}.zip
cp ncbi_dataset/data/${gcf}/* ${sp}
rm -r ncbi_dataset
rm README.md

## Deactivate env
conda deactivate
