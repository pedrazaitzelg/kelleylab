#!/bin/bash

#SBATCH --partition=256x44               # Partition/queue to run on (usually 128x24)
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

### for paralellizing each star run for SRA samples into a job array

### single star run for 1 pass of 1 sample of dwarf lemur

cd /hb/groups/kelley_training/itzel/data/transcriptomic

module load star

genome_dir=/hb/groups/kelley_training/itzel/anne/hibernation/star_out/dwarf_lemur #location of genome
fna=/hb/groups/kelley_training/itzel/anne/hibernation/data/genomic/dwarf_lemur/GCF_000165445.2_Mmur_3.0_genomic.fna     #location of fna file
trimmed_dir=/hb/groups/kelley_training/itzel/anne/hibernation/trimgalore_out/dwarf_lemur/white_adipose/trimgalore    #location of trimmed directory
mkdir -p dwarf_lemur/SRR5993015
cd dwarf_lemur/SRR599301

# Map Reads
if [ -f SRR599301_Log.final.out ]; then
    echo "already finished"
else
    STAR --genomeDir ${genome_dir} --runThreadN 4 --outFilterMultimapNmax 1 --twopassMode Basic --sjdbGTFfile /hb/groups/kelley_training/itzel/anne/hibernation/data/genomic/dwarf_lemur//genomic.gff \
        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./SRR599301_ \
        --readFilesIn ${trimmed_dir}/SRR5993015_pass_1_val_1.fq.gz  
fi
