#!/bin/bash
#SBATCH --partition=128x24
#SBATCH --time=01:00:00
#SBATCH --job-name=bbmerge_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user=igpedraz@ucsc.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --output=bbmerge_run.out
#SBATCH --error=bbmerge_run.err
#SBATCH --no-requeue


# module to run bbmerge located under bbtools
module load hb hb-gnu bbtools/bbtools-39.01

# this line will run the program for the 2 reads
# in1 = and in2 = set to locations of reads
bbmerge.sh in1=/hb/groups/kelley_training/itzel/data/anne/data/transcriptomic/dwarf_lemur_crossleyi/white_adipose/SRR5993015_pass_1.fastq.gz in2=/hb/groups/kelley_training/itzel/data/anne/data/transcriptomic/dwarf_lemur_crossleyi/white_adipose/SRR5993015_pass_2.fastq.gz out=SRR5993015_merged.fq
