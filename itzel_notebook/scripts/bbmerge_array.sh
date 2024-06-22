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
bbmerge.sh in1=read1.fq in2=read2.fq out=merged.fq
