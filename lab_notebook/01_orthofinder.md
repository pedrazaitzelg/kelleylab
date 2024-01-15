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
| | Run 1 (12 species) | Run 2 (11 species) | Run 3 (10 species, FINAL) |
| :---------------- | ------: | ----: | ----: |
| parallel blast step | 4h 32m | 1h 57m | 3h 14m |
| resumed orthofinder step (24 cores) | 1-19:53:54 | 1-11:42:11 | 1-10:46:52 |
| total time | ~2d | ~1d 14h | ~1d 14h |

- Run 1 (12 species)
- Run 2 (11 species, without house mouse)
- Run 3 (10 species, without little brown bat)

## Analysis
### Orthofinder stats
| | Run 1 (12 species) | Run 2 (11 species) | Run 3 (10 species, FINAL) |
| :---------------- | ------: | ----: | ----: |
| Number of species	| 12	| 11 | 10 |
| Number of genes | 613377	| 517185 | 474079 |
| Number of genes in orthogroups | 600010 | 505902 | 463078 |
| Number of unassigned genes | 13367 | 11283 | 11001 |
| Percentage of genes in orthogroups | 97.8 | 97.8 | 97.7 |
| Percentage of unassigned genes | 2.2 | 2.2 | 2.3 |
| Number of orthogroups | 28528 | 26689 | 26314 |
| Number of species-specific orthogroups | 4353 | 3249 | 3163 |
| Number of genes in species-specific orthogroups | 19057 | 13249 | 13102 |
| Percentage of genes in species-specific orthogroups | 3.1 | 2.6 | 2.8 |
| Mean orthogroup size | 21 | 19 | 17.6 |
| Median orthogroup size | 15 | 14 | 13 |
| G50 (assigned genes) | 33 | 28 | 26 |
| G50 (all genes) | 32 | 27 | 25 |
| O50 (assigned genes) | 5254 | 5280 | 5301 |
| O50 (all genes) | 5457 | 5482 | 5520 |
| Number of orthogroups with all species present | 13597 | 13660 | 14256 |
| Number of single-copy orthogroups | 646 | 877 | 1076 |

### Making a distribution of orthogroup size for OGs with all species present
#### Initial OrthoFinder run with 12 species:

The Orthogroups.GeneCount.tsv file in the Orthogroups results directory contains a table of per-species gene counts in each OG.
- path: /hb/groups/kelley_lab/anne/hibernation/orthofinder_run/orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out/Orthogroups/Orthogroups.GeneCount.tsv
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
- R Notebook for making distributions: [og_size_hist.Rmd](https://github.com/aanakamo/kelleylab_rotation/blob/main/scripts/og_size_hist.Rmd)

![Orthogroup size distribution](supporting_images/og_size_dist_reg_zoomed.png)

#### Redo this for OrthoFinder run with 11 species (no house mouse):

![Orthogroup size distribution](supporting_images/og_size_dist_reg_zoomed_11.png)

#### Redo this for final OrthoFinder run with 10 species (no little brown bat):

![Orthogroup size distribution](supporting_images/og_size_dist_reg_zoomed_10.png)

### Species trees
Made using [iTOL](https://itol.embl.de/upload.cgi)

11 species
~~~
(monito_del_monte:0.121445,(((brown_bear:0.00218596,black_bear:0.00243233)N4:0.0540776,(greater_horseshoe_bat:0.0526106,(brandts_bat:0.0169209,little_brown_bat:0.0151492)N7:0.0730869)N5:0.0161152)N2:0.0132531,(dwarf_lemur:0.0661784,((thirteen-lined_ground_squirrel:0.0114658,arctic_ground_squirrel:0.00594189)N8:0.0590534,(djungarian_hamster:0.0272538,syrian_hamster:0.0274553)N9:0.0848363)N6:0.0177884)N3:0.00985084)N1:0.121445)N0;
~~~
![Hibernating species tree](supporting_images/hibernating_species_tree.jpg)

10 species
~~~
(monito_del_monte:0.109891,((dwarf_lemur:0.0573769,((arctic_ground_squirrel:0.00552107,13_lined_ground_squirrel:0.00899323)N7:0.0511515,(djungarian_hamster:0.0246411,syrian_hamster:0.0260814)N8:0.0784609)N4:0.0163031)N2:0.00945774,((brandts_bat:0.0788643,greater_horseshoe_bat:0.0461365)N5:0.0139326,(black_bear:0.00225688,brown_bear:0.00222661)N6:0.0487507)N3:0.011697)N1:0.109891)N0;
~~~
![Hibernating species tree](supporting_images/hibernating_species_tree_10.jpg)
