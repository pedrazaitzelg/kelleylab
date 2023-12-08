#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=star                  # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=aanakamo@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=1                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-187]                  # array job

### for paralellizing each star run for SRA samples into a job array

cd /hb/groups/kelley_lab/anne/hibernation/star_out

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ../data/transcriptomic/species_tissue_sra_state.txt)
species=$(echo ${LINE} | awk '{ print $1; }')
tissue=$(echo ${LINE} | awk '{ print $2; }')
sra_acc=$(echo ${LINE} | awk '{ print $3; }')
state=$(echo ${LINE} | awk '{ print $4; }')

echo "running STAR for sra sample: ${sra_acc} (${species}, ${tissue}, ${state})"

genome_dir=../data/genomic/${species}
trimmed_dir=../trimgalore_out/${species}/${tissue}/trimgalore
mkdir -p ${species}/${tissue}

# Index genome for use with STAR
STAR --runMode genomeGenerate --genomeDir ${genome_dir} --genomeFastaFiles GCF_*_genomic.fna --sjdbGTFfile genomic.gff

# Map Reads
STAR --genomeDir ${genome_dir} --outFilterMultimapNmax 1 --twopassMode Basic --sjdbGTFfile genomic.gff \
        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./${species}/${tissue}/${sra_acc}_ \
        --readFilesIn ${trimmed_dir}/${sra_acc}_pass_1_val_1.fq ${trimmed_dir}/${sra_acc}_pass_2_val_2.fq
