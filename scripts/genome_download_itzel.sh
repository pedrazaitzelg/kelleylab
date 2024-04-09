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

cd /hb/groups/kelley_training/itzel/data/non_hibernators_genome

module load miniconda3.9

conda activate ncbi_datasets

while read line; do
    sp=$(echo ${line} | awk '{ print $1; }')
    #edited with -v RS+'\r\n' because of issues related to the script grabbing gcf from txt
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
