#!/bin/bash

## Login to the data transfer node (@hbfeeder.ucsc.edu)
## Run this script to download the genomic data for each hibernating species

cd /hb/groups/kelley_lab/anne/hibernation/data/genomic

source activate /hb/home/aanakamo/.conda/envs/ncbi_datasets

while read line; do
    sp=$(echo ${line} | awk '{ print $1; }')
    gcf=$(echo ${line} | awk '{ print $2; }')
    echo "*** downloading genomic data for ${sp}, accession ${gcf} ***"
    datasets download genome accession ${gcf} --include genome,protein,gff3,rna,cds --filename ${sp}.zip
done < species_gcf.txt

unzip *.zip

conda deactivate
