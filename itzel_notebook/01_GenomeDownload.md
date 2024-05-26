# Downloading Genomes from NCBI Database

https://www.ncbi.nlm.nih.gov/

1. Search ```Taxonomy``` to locate species of interest.
2. Select species and locate ```Genome```
3. From there, select from list of genomes available for that species.

```GCA_###``` refers to GenBank

```GCF_###``` refers to RefSeq 

"The Reference Sequence (RefSeq) collection provides a comprehensive, integrated, non-redundant, well-annotated set of sequences, including genomic DNA, transcripts, and proteins."

5. Download genome from NCBI using ```itzel_notebook/scripts/genome_download.sh``` and a corresponding ```species_gcf.txt``` file.


# Example of ```species_gcf.txt``` file
ringed_tailed_lemur  GCF_020740605.2

grey_mouse_lemur  GCF_000165445.2
