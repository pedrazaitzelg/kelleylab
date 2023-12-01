#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=5-00:00:00                # Max time for job to run
#SBATCH --job-name=orthofinder           # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=aanakamo@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=24               # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=120G                       # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

# cd /hb/groups/kelley_lab/anne/hibernation
# mkdir -p orthofinder_run
# cd orthofinder_run

# mkdir -p orthofinder_in    ## directory to contain the input faa files
# cat ../data/genomic/species_gcf.txt| awk '{ print $1; }' | while read sp; do
#     cp ../data/genomic/${sp}/protein.faa orthofinder_in/${sp}.faa
# done

### process the faa files to include species names in each gene
# cd orthofinder_in
# source activate /hb/home/aanakamo/.conda/envs/biopython
# python ~/kelley_lab_rotation/scripts/prep_faa_for_orthofinder.py
# conda deactivate
# rm !(*.prepped.faa)
# cd ..

cd /hb/groups/kelley_lab/anne/hibernation/orthofinder_run
source activate /hb/home/aanakamo/.conda/envs/orthofinder
### Testing orthofinder on the provided example data
#orthofinder -f ExampleData -t 24 -a 5 -M msa -A mafft -T fasttree -o ExampleData_out -S diamond_ultra_sens
### One-go orthofinder command
#orthofinder -f orthofinder_in -t 24 -a 5 -M msa -A mafft -T fasttree -o orthofinder_out -S diamond_ultra_sens

### run blast step separately
orthofinder -op -S diamond_ultra_sens -f orthofinder_in -n out -o orthofinder_out | grep "diamond blastp" > jobqueue
mv jobqueue jobqueue_old
shuf jobqueue_old > jobqueue
sbatch ~/kelleylab_rotation/scripts/orthofinder_blast_array.sh

### after separate blast step
orthofinder -M msa -A mafft -T fasttree -t 24 -a 5 -n out -b orthofinder_out/Results_out/WorkingDirectory

conda deactivate
