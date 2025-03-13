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

#### WORKS !!! BUT SET TO ONLY ONE LINE SO NOW 
#### NEED TO GENERATE ARRAY TO GO THROUGH LIST OF VALUES


## Filter .frq files generated from vcf_tools to remove
## things fixed in both populations an doutput only the differences

#set working directory
cd /hb/groups/kelley_training/itzel/population_bears_proj24/population_stats/allele_freq/temp2

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

chrom=$( awk ' NR==2{print $1} ' chrom.txt)  #creates variable set to a chrom name

# search file for chromosome and output new file with only that chromosome position
grep ${chrom} ${allele_file1} > ${chrom}_${location_1}.txt
grep ${chrom} ${allele_file2} > ${chrom}_${location_2}.txt

## sets variable to values in second column
pos1=$(awk 'NR==1 {print $2}' ${chrom}_${location_1}.txt)
freq1=$(awk 'NR==1 {print $5}' ${chrom}_${location_1}.txt)

pos2=$(awk 'NR==1 {print $2}' ${chrom}_${location_2}.txt)
freq2=$(awk 'NR==1 {print $5}' ${chrom}_${location_2}.txt)

## output difference in freq between both files
diff=$(awk -v f1="$freq1" -v f2="$freq2" 'BEGIN {print f1 - f2}')  # -v sets variable to preset

echo 'CHROM = ' ${chrom}
echo 'Position =' ${pos1}
echo 'Location 1:' ${location_1}
echo 'Location 2:' ${location_2}
echo 'Difference:' ${freq1} '-' ${freq2} '=' ${diff}
#echo ${chrom}  ${pos1}  ${diff}

#remove all temp files
#rm ${chrom}_1.txt
#rm ${chrom}_2.txt
