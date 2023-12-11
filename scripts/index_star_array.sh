#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=2-00:00:00                # Max time for job to run
#SBATCH --job-name=star                  # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=aanakamo@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=24                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=30G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[5]                   # array job

### for paralellizing each star run for SRA samples into a job array

#cd /hb/groups/kelley_lab/anne/hibernation/star_out
cd /hb/scratch/aanakamo/kelleylab_rotation/star_tmp

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /hb/scratch/aanakamo/kelleylab_rotation/star_tmp/species_gcf.txt)
species=$(echo ${LINE} | awk '{ print $1; }')

echo "running STAR indexing for: ${species}"

genome_dir=/hb/groups/kelley_lab/anne/hibernation/data/genomic/${species}
fna=$(basename ${genome_dir}/GCF_*_genomic.fna)
#trimmed_dir=../trimgalore_out/${species}/${tissue}/trimgalore
mkdir -p ${species}
cd ${species}

# Index genome for use with STAR
#STAR --runMode genomeGenerate --runThreadN 24 --genomeDir . --genomeFastaFiles ${genome_dir}/${fna} --sjdbGTFfile ${genome_dir}/genomic.gff
STAR --runMode genomeGenerate --runThreadN 24 --genomeDir . --genomeFastaFiles ${genome_dir}/${fna} --sjdbGTFfile ${genome_dir}/genomic.gff --limitGenomeGenerateRAM 123560700863
