#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=2-00:00:00                # Max time for job to run
#SBATCH --job-name=featCounts            # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=aanakamo@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-18]                   # array job (18, number of matrices)

### for paralellizing each featurecounts run for count matrices into a job array

cd /hb/scratch/aanakamo/kelleylab_rotation/featurecounts_tmp

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p matrices_list.txt)
species=$(echo ${LINE} | awk '{ print $1; }')
tissue=$(echo ${LINE} | awk '{ print $2; }')

mkdir -p ${species}/${tissue}
cd ${species}/${tissue}

echo "running featureCounts for count matrix: ${species}, ${tissue}"
bam_files=/hb/scratch/aanakamo/kelleylab_rotation/star_tmp/${species}/${tissue}/*/*.sortedByCoord.out.bam
saf_file=../${species}.OG.allGeneInfo.saf

# Quantify gene-level counts using featureCounts
source activate /hb/home/aanakamo/.conda/envs/featurecounts
featureCounts -p -F 'SAF' -T 8 -t exon -g gene_id -a ${saf_file} -o ./${species}.${tissue}.featurecounts ${bam_files}
conda deactivate
