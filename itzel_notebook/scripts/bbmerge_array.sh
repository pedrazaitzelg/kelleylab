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
#SBATCH --array=[3-14]


# module to run bbmerge located under bbtools
module load hb hb-gnu bbtools/bbtools-39.01

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p species_tissue_sra_state.txt)
species=$(echo ${LINE} | awk '{ print $1; }')
tissue=$(echo ${LINE} | awk '{ print $2; }')
sra_acc=$(echo ${LINE} | awk '{ print $3; }')
state=$(echo ${LINE} | awk '{ print $4; }')

sra_dir=/hb/groups/kelley_training/itzel/data/anne/data/transcriptomic/${species}/${tissue}

echo "running bbmerge for sra sample: ${sra_acc} (${species}, ${tissue}, ${state})"

mkdir -p ${species}/${tissue}
cp bbmerge.sh ${species}/${tissue}/
cd  ${species}/${tissue}

# this line will run the program for the 2 reads
# in1 = and in2 = set to locations of reads
../../bbmerge.sh in1=${sra_dir}/${sra_acc}_pass_1.fastq.gz in2=${sra_dir}/${sra_acc}_pass_2.fastq.gz out=${sra_acc}_merged.fq
