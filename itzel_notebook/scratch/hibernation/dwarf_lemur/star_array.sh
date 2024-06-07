#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on (usually 128x24)
#SBATCH --time=2-00:00:00                # Max time for job to run
#SBATCH --job-name=star                  # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --output=star.out
#SBATCH --error=star.err
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=4                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=40G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-22]                 # array job

### for paralellizing each star run for SRA samples into a job array

cd /hb/groups/kelley_training/itzel/data/transcriptomic/star_out

module load star

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /hb/groups/kelley_training/itzel/data/transcriptomic/species_tissue_sra_state.txt)
species=$(echo ${LINE} | awk '{ print $1; }')
tissue=$(echo ${LINE} | awk '{ print $2; }')
sra_acc=$(echo ${LINE} | awk '{ print $3; }')
state=$(echo ${LINE} | awk '{ print $4; }')

echo "running STAR for sra sample: ${sra_acc} (${species}, ${tissue}, ${state})"

genome_dir=/hb/groups/kelley_training/itzel/data/genomic/lesser_dwarf_lemur #location of genome
fna=/hb/groups/kelley_training/itzel/data/genomic/lesser_dwarf_lemur/GCF_*_genomic.fna
trimmed_dir=/hb/groups/kelley_training/itzel/anne/hibernation/trimgalore_out/dwarf_lemur/white_adipose/trimgalore
mkdir -p ${species}/${tissue}/${sra_acc}
cd ${species}/${tissue}/${sra_acc}

# Map Reads
if [ -f ${sra_acc}_Log.final.out ]; then
    echo "already finished"
else
 STAR --genomeDir ${genome_dir} --runThreadN 20 --outFilterMultimapNmax 1 --twopassMode Basic --sjdbGTFfile ${genome_dir}/genomic.gff \
        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./${sra_acc}_ \
        --readFilesIn ${trimmed_dir}/${sra_acc}_pass_1_val_1.fq.gz ${trimmed_dir}/${sra_acc}_pass_2_val_2.fq.gz
fi
