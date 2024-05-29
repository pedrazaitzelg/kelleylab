#!/bin/bash
#SBATCH --job-name=genomic_data     #name of job
#SBATCH --partition=128x24          #partition for job to run
#SBATCH --mail-type=ALL              #updates from hb
#SBATCH --mail-user=igpedraz@ucsc.edu  
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output=genomic_data.out    #output files
#SBATCH --error=genomic_data.err    #error files
#SBATCH --mem=250M

## Run this script to download the genomic data for each species using species_gcf.txt

cd /hb/groups/kelley_training/itzel/new_data/genomic

module load miniconda3.9

conda activate ncbi_datasets

while read line; do
    sp=$(echo ${line} | awk '{ print $1; }')
    gcf=$(echo ${line} | awk -v RS='\r\n' '{ print $2; }')
    echo "*** downloading genomic data for ${sp}, accession ${gcf} ***"
    datasets download genome accession "${gcf}" --include genome,protein,gff3,rna,cds --filename ${sp}.zip

    # If the previous command does not succeed, break out of the while loop immediately.
    if [ $? -ne 0 ]; then break; fi

    unzip ${sp}.zip
    mkdir -p ${sp}
    cp ncbi_dataset/data/${gcf}/* ${sp}
    rm -r ncbi_dataset       
    rm README.md
done < species_gcf.txt

rm *.zip

conda deactivate
