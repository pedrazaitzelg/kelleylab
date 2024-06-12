#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=fastqc                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=bowtie2.out            # Standard output and error log
#SBATCH --error=bowtie2.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[3]                  # array job

cd /hb/groups/kelley_training/itzel/data/fastq_screen/index_bowtie

module load bowtie/bowtie2-2.3.2

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /hb/groups/kelley_training/itzel/data/fastq_screen/index_bowtie/species_gcf.txt)
species=$(echo ${LINE} | awk '{ print $1; }')
gcf=$(echo ${LINE} | awk '{ print $2; }')

genome_dir=/hb/groups/kelley_training/itzel/data/fastq_screen/index_bowtie/genomes/${species}  #location of genomic data
fna=$(basename ${genome_dir}/GCF_*_genomic.fna)    #location of .fna files  |   typically GCF files bute in not RefSeq then edit to GCA

echo "running bowtie index for: ${species} ${gcf}"

mkdir -p ${species}
cd ${species}

bowtie2-build genomes/${species}/${fna} ${species}    #location of fna file followed by index name for species
