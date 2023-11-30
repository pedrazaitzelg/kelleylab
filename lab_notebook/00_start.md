# Notes for getting started

[Markdown cheatsheet](https://www.markdownguide.org/cheat-sheet/)

Configure remote git repo (in home directory /hb/home/aanakamo): 

`git clone https://github.com/aanakamo/kelleylab_rotation.git`

`git config user.name "aanakamo"` 

`git config user.email "aanakamo@ucsc.edu"` 

`git remote set-url origin https://{personal-access-token}@github.com/aanakamo/kelleylab_rotation.git` 

Slurm on [Hummingbird](https://hummingbird.ucsc.edu/) 

`ssh aanakamo@hb.ucsc.edu`

`cd /hb/groups/kelley_lab/anne`

- running jobs: [tutorial](https://hummingbird.ucsc.edu/documentation/creating-scripts-to-run-jobs/)
- data transfer: `sftp aanakamo@hbfeeder.ucsc.edu`
- interactive job:
    `salloc --partition=128x24 --time=02:00:00 --mem=10G --ntasks=1 --cpus-per-task=1`

    `export | grep SLURM`

    `ssh $SLURM_NODELIST`

    `# run stuff`

    `exit`

    `exit`

Loading modules:

`module load fastqc`

`module load sratoolkit`

Installing software with conda:

`conda install -c bioconda orthofinder`

`conda install -c conda-forge ncbi-datasets-cli`

