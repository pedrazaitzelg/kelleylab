#!/bin/bash

## Login to the data transfer node (@hbfeeder.ucsc.edu)
## no longer need to enter transfer mode - can run on regular @hb.ucsc.edu
## Run this script to download the genomic data for each hibernating species

cd /hb/groups/kelley_training/itzel/data/genomic

source activate /hb/home/igpedraz/.conda/envs/ncbi_datasets

while read line; do
    sp=$(echo ${line} | awk '{ print $1; }')
    gcf=$(echo ${line} | awk '{ print $2; }')
    echo "*** downloading genomic data for ${sp}, accession ${gcf} ***"
    datasets download genome accession ${gcf} --include genome,protein,gff3,rna,cds --filename ${sp}.zip
    unzip ${sp}.zip
    mkdir -p ${sp}
    cp ncbi_dataset/data/${gcf}/* ${sp}
    rm -r ncbi_dataset
    rm README.md
done < species_gcf.txt

rm *.zip

conda deactivate
