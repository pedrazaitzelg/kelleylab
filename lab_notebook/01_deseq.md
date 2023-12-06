# Deseq2 Notes

[Deseq documentation](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) 

[Blair's RNA-Seq analysis github tutorial](https://github.com/blairperry/midhib_feeding_uarctos#1-quality-trimming-mapping-and-processing-of-rna-seq-data)

## Downloading transcriptomic data

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
- script: [download_transcrptomic_array.sh]()

