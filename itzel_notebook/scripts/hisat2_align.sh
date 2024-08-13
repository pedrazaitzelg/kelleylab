#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=hisat_run                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=hisat_run.out            # Standard output and error log
#SBATCH --error=hisat_run.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

#indexing reference genome for hisat2
cd /hb/groups/kelley_training/itzel/data/hisat2/test_run

#load module
module load hisat/2.1.0

#run hisat alignment
# hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]
# -q indicates using fastq files
hisat2 -q -x ../index_genomes/dwarf_lemur/ -1 ../transcriptomic_data/dwarf_lemur_crossleyi/white_adipose/SRR5993015_pass_1.fastq.gz  -2 ../transcriptomic_data/dwarf_lemur_crossleyi/white_adipose/SRR5993015_pass_2.fastq.gz -S SRR5993015.sam
