# Navigating Hummingbird

ssh igpedraz@hb.ucsc.edu

**home directory** /hb/groups/kelley_training/itzel

# NCBI_Datasets

Install using conda
The NCBI Datasets CLI tools are available as a conda package that includes both datasets and dataformat.

First, create a conda environment: ```conda create -n ncbi_datasets```

Then, activate your new environment: ```conda activate ncbi_datasets```

Finally, install the datasets conda package: ```conda install -c conda-forge ncbi-datasets-cli```

# Modules

```module load star``` to map data onto genomes using STAR

```module load condaminiconda3.9``` when downloading genome data from NCBI

```module load fastqc``` for creating FastQC files from downloaded genomic data

```module load trimgalore``` for trimming genomic data from fastqc files

# Useful Commands

```pwd``` print working directory

```mkdir``` make new directory

```nano``` create/edit file

```cd``` move between directories

```ls``` list what is in the current directory

```ls -alt``` lists hidden files, last data modified, in alphabetical order
