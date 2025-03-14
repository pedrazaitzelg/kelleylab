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
#SBATCH --array=[1-1]


## Filter .frq files generated from vcf_tools to remove
## things fixed in both populations an doutput only the differences

#set working directory
cd /hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/allele_freq/temp

##pairwise comparisons
#set variables
# Set directory
dir=/hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/indv_stats

# Read locations from the file into an array using mapfile function
# acess each line by echo "${locations[0]}"
mapfile -t locations < locations.txt

# Loop through each pair of locations
for ((i=0; i<${#locations[@]}; i++)); do  #sets 'i' as first location in array; '++' increases variable by 1 after each iteration
    for ((j=i+1; j<${#locations[@]}; j++)); do    #sets j as second location in array after the first
        location_1=${locations[i]}
        location_2=${locations[j]}

        echo "Comparing $location_1 to $location_2..."

        # Define allele frequency file paths
        allele_file1="${dir}/${location_1}/${location_1}.frq"
        allele_file2="${dir}/${location_2}/${location_2}.frq"

        # Define output file (tab-delimited format)
        output_file="${location_1}_${location_2}.txt"
        echo -e "CHROM"\t"POS"\t"Diff" > "$output_file"  # Adds header

        # Extract chromosome list 
        chrom_list=$(awk '{print $1}' "$allele_file1" | sort | uniq)

        # Loop through each chromosome
        for chr in $chrom_list; do
            if [ "$chr" != "CHROM" ]; then
                # Extract position and frequency information
                pos1_list=$(awk -v chr="$chr" '$1 == chr {print $2}' "$allele_file1")
                pos2_list=$(awk -v chr="$chr" '$1 == chr {print $2}' "$allele_file2")

                # Loop through positions in location_1
                for p1 in $pos1_list; do
                    # Check if the position exists in location_2
                    if echo "$pos2_list" | grep -wq "$p1"; then
                        # Extract corresponding allele frequencies
                        freq1=$(awk -v chr="$chr" -v p="$p1" '$1 == chr && $2 == p {print $5}' "$allele_file1")
                        freq2=$(awk -v chr="$chr" -v p="$p1" '$1 == chr && $2 == p {print $5}' "$allele_file2")

                        # Compute frequency difference
                        diff=$(awk -v f1="$freq1" -v f2="$freq2" 'BEGIN {print f1 - f2}')

                        # Append results to output file
                        echo -e "$chr\t$p1\t$diff" >> "$output_file"
                    fi
                done
            fi
        done

        echo "Results saved to $output_file (tab-delimited)"
    done
done



