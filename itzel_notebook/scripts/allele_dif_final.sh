#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=allele_filter                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=igpedraz@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-13]


## Filter .frq files generated from vcf_tools to remove
## things fixed in both populations an doutput only the differences

#set working directory
cd /hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/allele_freq/temp

##pairwise comparisons
## start with ABC and Alaska

#set variables
#LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p locations.txt)
dir=/hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/indv_stats

## sets location name from locations.txt
location_1=$( awk ' NR==1{print $1} ' locations.txt)
location_2=$( awk ' NR==2{print $1} ' locations.txt)


allele_file1=${dir}/${location_1}/${location_1}.frq  #location of files
allele_file2=${dir}/${location_2}/${location_2}.frq  #location of files

awk ' {print $1; } ' ${allele_file1} | sort | uniq > chrom.txt  #generate txt file with list of names

# sets chrom to equal list of chrom names to run through
# outputs all chromosome files for each location
chrom=$( awk ' {print $1;} ' chrom.txt)

for chr in ${chrom}; do
        if [ "$chr" != "CHROM" ]; then
                # Search for chromosome and output new file with only that chromosome's positions
                grep "$chr" "$allele_file1" > "${chr}_${location_1}.txt"
                grep "$chr" "$allele_file2" > "${chr}_${location_2}.txt"
        fi
done


# now go through each position
# Set output file dynamically based on location names
output_file="${location_1}_${location_2}.txt"
echo -e "Chromosome\tPosition\tDiff" > "$output_file"  # Add header with tabs

# for each chromosome, go through each position 
for chr in ${chrom}; do
    if [ "$chr" != "CHROM" ]; then
        # Set position and frequency variables
        pos1_list=$(awk '{print $2}' ${chr}_${location_1}.txt)
        pos2_list=$(awk '{print $2}' ${chr}_${location_2}.txt)

        # Loop over positions in location_1
        for p1 in ${pos1_list}; do
            # Check if the position exists in location_2
            if echo "${pos2_list}" | grep -wq "${p1}"; then
                # Extract corresponding frequencies for that position
                freq1=$(awk -v p="$p1" '$2 == p {print $5}' ${chr}_${location_1}.txt)
                freq2=$(awk -v p="$p1" '$2 == p {print $5}' ${chr}_${location_2}.txt)

                # Compute frequency difference
                diff=$(awk -v f1="$freq1" -v f2="$freq2" 'BEGIN {print f1 - f2}')

                # set output file
                output_file="${location_1}_${location_2}.txt"
                # Append result to output file
                echo -e "$chr\t$p1\t$diff" >> "$output_file"
            fi
        done
    fi
done



              
               
