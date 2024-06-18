#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=getSRA                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=1                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=SRArun_%j.out            # Standard output and error log
#SBATCH --error=SRArun_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-11]                  # array job

### for paralellizing each SRA sample download into a job array

#directory location for downloading to take place
cd /hb/groups/kelley_training/itzel/data/transcriptomic 

#activates module environment
module load  miniconda3.9

#activates module but may not be necessary; included anyway
conda activate ncbi_datasets

#module needed to use prefetch and fastq-dump commands
module load  hb  hb-gnu  sratoolkit/sratoolkit-3.0.10

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p species_tissue_sra_state.txt)
species=$(echo ${LINE} | awk '{ print $1; }')
tissue=$(echo ${LINE} | awk '{ print $2; }')
sra_acc=$(echo ${LINE} | awk '{ print $3; }')
state=$(echo ${LINE} | awk '{ print $4; }')

echo "downloading sra sample: ${sra_acc} (${species}, ${tissue}, ${state})"

prefetch --output-directory ${species}/${tissue} ${sra_acc}        # for more space --max-size followed by #G
fastq-dump --outdir ${species}/${tissue} --skip-technical --readids --read-filter pass \
    --gzip --dumpbase --split-3 --clip ${species}/${tissue}/${sra_acc}/${sra_acc}.*
rm -r ${species}/${tissue}/${sra_acc}

conda deactivate
