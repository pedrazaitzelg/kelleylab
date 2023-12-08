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
Available on Hummingbird as a module (is a conda env): `module load trimgalore`

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
    - Blair: this was based on the data in the particular study; choose something and stick with it for consistency, don't need to trim samples differently
    - Looking at my fastqc reports, all samples (except those from just one study), could benefit from trimming the first 8b at the 5' end. trim_galore already handles low quality bases at the 3' end, and removes reads with overall low quality, so this takes care of the other cases. I'll use `--clip_R1 8` and `--clip_R2 8` for now and look at the fastqc reports after trimming 

my command: `trim_galore --paired -q 20 --fastqc --fastqc_args "--nogroup --outdir [fastqc/out/dir]" --stringency 5 --illumina --length 50 -o [trimgalore/out/dir] --clip_R1 8 --clip_R2 8 [path/to/read1] [path/to/read2]`
- array job script: [trimgalore_array.sh](https://github.com/aanakamo/kelleylab_rotation/blob/main/scripts/trimgalore_array.sh)

### Mapping reads with STAR
Available on Hummingbird as a module: `module load star`

Blair's commands:
~~~
# Index genome for use with STAR
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ./star_reference --genomeFastaFiles GCF_023065955.1_UrsArc1.0_genomic.fna --sjdbGTFfile GCF_023065955.1_UrsArc1.0_genomic.gff

# Map Reads
STAR --genomeDir ./star_reference/ --runThreadN 8 --outFilterMultimapNmax 1 --twopassMode Basic --sjdbGTFfile GCF_023065955.1_UrsArc1.0_genomic.gff --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./star_mapped/[output_prefix] --readFilesIn [path/to/file1] [path/to/file2]
~~~
- `--outFilterMultimapNmax 1` - Default 10, "maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as "mapped to too many loci" in the Log.final.out." So, 1 is strict. **If there are multiple isoforms, would this cause issues?**
    - Seems like gff files don't have duplicate exons? Need to check on this further...
- `--twopassMode Basic` - "basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly"
- `--readFilesCommand zcat` - command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout. For example: zcat - to uncompress .gz files, bzcat - to uncompress .bz2 files, etc." 
- `--outSAMtype BAM SortedByCoordinate` - output BAM file, "sorted by coordinate. This option will allocate extra memory for sorting which can be specified by --limitBAMsortRAM."

my commands:
~~~
# Index genome for use with STAR
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir [path/to/genomic] --genomeFastaFiles GCF_*_genomic.fna --sjdbGTFfile genomic.gff

# Map Reads
STAR --genomeDir [path/to/genomic] --runThreadN 8 --outFilterMultimapNmax 1 --twopassMode Basic --sjdbGTFfile genomic.gff --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix [star/output/dir/prefix] --readFilesIn [path/to/file1] [path/to/file2]
~~~
- array job script: [star_array.sh](https://github.com/aanakamo/kelleylab_rotation/blob/main/scripts/star_array.sh)

### Create count matrices with featureCounts
Blair's commands:
~~~
# Convert gff to gtf w/ AGAT
agat_convert_sp_gff2gtf.pl --gff GCF_023065955.1_UrsArc1.0_genomic.gff -o GCF_023065955.1_UrsArc1.0_genomic.gtf

# Convert gtf to SAF using custom python script 
python2 ./utility_scripts/gtf_to_saf.allGeneInfo.py GCF_023065955.1_UrsArc1.0_genomic.gtf exon
# Output file is titled GCF_023065955.1_UrsArc1.0_genomic.allGeneInfo.saf

# Quantify gene-level counts using featureCounts
featureCounts -p -F 'gtf' -T 8 -t exon -g gene_id -a GCF_023065955.1_UrsArc1.0_genomic.allGeneInfo.saf -o ./postDexExperiment_08.10.22.txt [path/to/star_mapped]/*.sortedByCoord.out.bam
~~~
- `-p` - "If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads; single-end reads are always counted as reads." This just indicates the reads are paired end.
- `-F 'gtf'` - "Specify format of the provided annotation file. Acceptableformats include 'GTF' (or compatible GFF format) and 'SAF'. 'GTF' by default. For SAF format, please refer to Users Guide."
- `-t exon` - "Specify feature type(s) in a GTF annotation ... 'exon' by default. Rows in the annotation with a matched feature will be extracted and used for read mapping. "
- `-g gene_id` - "Specify attribute type in GTF annotation. 'gene_id' by default. Meta-features used for read counting will be extracted from annotation using the provided value."

## The Deseq step

