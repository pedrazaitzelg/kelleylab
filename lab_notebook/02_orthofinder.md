# OrthoFinder Notes
OrthoFinder will be used to find orthologous genes across hibernating species, which will help to implement cross-species differential expression analysis.

## Installing OrthoFinder
`conda create -n orthofinder`

`conda activate orthofinder`

`conda install -c bioconda orthofinder`
- orthofinder conda package came with the wrong version of diamond, so manually installed diamond and replaced the executable in the conda env directory
~~~
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
cp diamond /hb/home/aanakamo/.conda/envs/orthofinder/bin/
~~~

## Running OrthoFinder
#### Scripts used and description:
- [download_genomic.sh](https://github.com/aanakamo/kelleylab_rotation/blob/main/scripts/download_genomic.sh) - downloads the genomic data for hibernating species of interest, including the .faa protein files needed for OrthoFinder
- [orthofinder_run.sh](https://github.com/aanakamo/kelleylab_rotation/blob/main/scripts/orthofinder_run.sh) - shows the overall process of running OrthoFinder on Hummingbird cluster:
    - sets up the necessary directories
    - uses [prep_faa_for_orthofinder.py](https://github.com/aanakamo/kelleylab_rotation/blob/main/scripts/prep_faa_for_orthofinder.py) to process each .faa file for orthofinder by adding the species name to each gene and removing stop codons
    - the blast step of orthofinder must be parallelized to avoid running out of memory on Hummingbird, so each command is obtained from the orthofinder log
    - then [orthofinder_blast_array.sh](https://github.com/aanakamo/kelleylab_rotation/blob/main/scripts/orthofinder_blast_array.sh) array job can be run separately to parallelize the blast commands
    - when the blast jobs are finished, orthofinder can be resumed from the blast results

#### Running time:
| | Run 1 (12 species) | Run 2 (11 species, without house mouse) |
| :---------------- | ------: | ----: |
| parallel blast step | 4h 32m | 1h 57m |
| resumed orthofinder step | 1-19:53:54 |  |
| total time | ~2 days |  |

## Analysis
### Orthofinder stats
| | Run 1 (12 species) | Run 2 (11 species, without house mouse) |
| :---------------- | ------: | ----: |
| Number of species	| 12	| 
| Number of genes | 613377	| 
| Number of genes in orthogroups | 600010 | 
| Number of unassigned genes | 13367 |
| Percentage of genes in orthogroups | 97.8 |
| Percentage of unassigned genes | 2.2 |
| Number of orthogroups | 28528 |
| Number of species-specific orthogroups | 4353 |
| Number of genes in species-specific orthogroups | 19057 |
| Percentage of genes in species-specific orthogroups | 3.1 |
| Mean orthogroup size | 21 |
| Median orthogroup size | 15 |
| G50 (assigned genes) | 33 |
| G50 (all genes) | 32 |
| O50 (assigned genes) | 5254 |
| O50 (all genes) | 5457 |
| Number of orthogroups with all species present | 13597 |
| Number of single-copy orthogroups | 646 |

### Making a distribution of orthogroup size for OGs with all species present
#### Initial OrthoFinder run with 12 species:

The Orthogroups.GeneCount.tsv file in the Orthogroups results directory contains a table of per-species gene counts in each OG.
- path: /hb/groups/kelley_lab/anne/hibernation/orthofinder_run/old_orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out/Orthogroups/Orthogroups.GeneCount.tsv
- want to condition on each value in a row being non-zero
- end up with a table with cols: OG, total_gene_count
- quick python script to do this:
~~~
import sys

for line in sys.stdin:
    if "OG0" in line:
        lst = line.strip().split()
        og = lst[0]; total = lst[-1]
        include = True
        for i in range(1,len(lst)):
            if int(lst[i]) == 0:
                include = False
        if include:
            print(og + "\t" + total)
~~~
![Orthogroup size distribution](og_size_dist_reg_zoomed.png)

#### Redo this for final OrthoFinder run with 11 species (no house mouse):
