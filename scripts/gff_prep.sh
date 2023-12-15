#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=prepGFF               # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=aanakamo@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=2                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[9]                   # array job (1-10)

### for paralellizing preparation of gff files for input to featureCounts into a job array

cd /hb/scratch/aanakamo/kelleylab_rotation/featurecounts_tmp

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /hb/groups/kelley_lab/anne/hibernation/data/genomic/species_gcf.txt)
species=$(echo ${LINE} | awk '{ print $1; }')

echo "prepping gff for: ${species}"

gff_file=/hb/groups/kelley_lab/anne/hibernation/data/genomic/${species}/genomic.gff
og_dir=/hb/groups/kelley_lab/anne/hibernation/orthofinder_run/orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out/Orthogroups
col_num=$(grep ${species} species_index.txt | awk '{ print $2; }')
mkdir -p ${species}
cd ${species}

### parse orthogroup tsv file, output has cols: "protein_id", "Orthogroup"
python2 ~/kelleylab_rotation/scripts/blairs_gff_scripts/parseOrthogroups.py ${og_dir}/Orthogroups.tsv ${species}.prepped ${col_num}
### make a mapping between "protein_id" and "GeneID"
python ~/kelleylab_rotation/scripts/og_proteinid.py ${gff_file} ${og_dir}/Orthogroups.${species}.prepped.tsv

### add "protein_id", "Orthogroup" info to the gff file
python2 ~/kelleylab_rotation/scripts/blairs_gff_scripts/annotateGFF.py ${gff_file} ${og_dir}/Orthogroups.${species}.prepped.geneID.tsv ${species}.OG.gff

### convert GFF to GTF w/ AGAT
module load agat
agat_convert_sp_gff2gtf.pl --gff ${species}.OG.gff -o ${species}.OG.gtf
module unload agat

### Convert GTF to SAF

