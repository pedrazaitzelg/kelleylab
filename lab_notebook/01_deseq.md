# DESeq2 Notes

Deseq2 will be used to find differentially expressed genes within and across hibernating species. Some resources to get started:

[Deseq documentation](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) 

[Blair's RNA-Seq analysis github tutorial](https://github.com/blairperry/midhib_feeding_uarctos#1-quality-trimming-mapping-and-processing-of-rna-seq-data)

## Pre-Deseq steps

### Downloading transcriptomic data

Following instructions from here: [Batch downloading FASTQ files using the SRA toolkit, fastq-dump, and Python](https://erilu.github.io/python-fastq-downloader/)
- logged onto hummingbird feeder
- First trying on SRR6131236, a brown bear (fat, active) sample
~~~
cd
prefetch SRR6131236
fastq-dump --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR6131236/SRR6131236.sra
rm -r SRR6131236
~~~
- now writing an array job to do this for all samples
- /hb/groups/kelley_lab/anne/hibernation/data/transcriptomic/species_tissue_sra_state.txt - contains a full list of the samples to download (cols: species, tissue, SRA accession, state)
- script: [download_transcrptomic_array.sh](https://github.com/aanakamo/kelleylab_rotation/blob/main/scripts/download_transcrptomic_array.sh)

### Assessing quality with fastqc
- First trying on SRR6131236 again:
~~~
fastqc --extract --outdir fastqc_out SRR6131236_pass_1.fastq SRR6131236_pass_2.fastq
~~~
- guide for interpreting fastqc html output: [Introduction to RNA-Seq using high-performance computing](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html)

array job script: [fastqc_array.sh](https://github.com/aanakamo/kelleylab_rotation/blob/main/scripts/fastqc_array.sh)

### Trimming with trim_galore
Available on Hummingbird as a module (is a conda env)
`module load trimgalore`

Blair's command:
`trim_galore --paired -q 20 --fastqc --fastqc_args "--noextract --nogroup --outdir 2_TrimGalore/fastqc/" --stringency 5 --illumina --length 50 -o trimmed_reads/ --clip_R1 12 --clip_R2 12 [path/to/read1] [path/to/read2]`
- `--paired` - specify the files are paired-end
- `-q 20` - "Trim low-quality ends from reads in addition to adapter removal ... Default Phred score: 20"
- `--fastqc` - run fastqc after trimming
- `--fastqc_args "--noextract --nogroup --outdir 2_TrimGalore/fastqc/"` - specify arguments for running fastqc
- `--stringency 5` - "Overlap with adapter sequence required to trim a sequence. Defaults to a very stringent setting of 1, i.e. even a single bp of overlapping sequence will be trimmed off from the 3' end of any read." So, 5 would be less stringent
- `--illumina` - "Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter 'AGATCGGAAGAGC' instead of the default auto-detection of adapter sequence." (I checked that all the samples are illumina)
- `--length 50` - "Discard reads that became shorter than length INT because of either quality or adapter trimming. A value of '0' effectively disables this behaviour. Default: 20 bp." So, 50 is a bit more strict
- `--clip_R1 12` and `--clip_R2 12` - "Instructs Trim Galore to remove *int* bp from the 5' end of ... This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at the 5' end." **Do I need this?**
