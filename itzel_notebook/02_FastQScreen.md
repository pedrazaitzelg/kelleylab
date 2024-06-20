# FastQ Screen to look for potential contamination in dataset

# 01 - Downloading fastq screen
`wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.tar.gz`
`tar -xzf v0.15.3.tar.gz`

#### downloading Bowtie2 and adding bowtie to path
`wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip/download`
`unzip download`
`export PATH=$HOME/tools/bowtie2/bowtie2-2.4.2-sra-linux-x86_64:$PATH`

#### adding bismark to path, for bisulfite data analysis only (example from Tina)
`export PATH=/hb/groups/kelley_lab/tina/mytilus/Bismark-master/bismark:$PATH`


# 02 - Index Reference Genomes with Bowtie2
`module load bowtie/bowtie2-2.3.2`
`bowtie2-build ../genomes/${species}/GCF_* ${species}`

# 03 - Running fastq screen

#### (1) download the test dataset and run test
`wget https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_test_dataset.tar.gz`
`tar -xzf fastq_screen_test_dataset.tar.gz`

#### sample: /hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim/concat_lanes/12O-F_S66_R1_merged.fq.gz

#### before running fastq screen, need a configuration file (indicating path of genomes) in the directory 
`cd /hb/groups/kelley_training/itzel/fastqc_screen/FastQ-Screen-0.15.3`
`cat fastq_screen.conf (see fastq_screen.config)`

#### running Fastq screen for bisulfite data (example from Tina)
`FastQ-Screen-0.15.3/fastq_screen --bisulfite  /hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim/concat_lanes/12O-F_S66_R1_merged.fq.gz`
