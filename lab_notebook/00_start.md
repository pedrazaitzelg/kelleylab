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

### Loading modules:

`module load fastqc`

`module load sratoolkit`

### Installing software with conda:

`conda install -c bioconda orthofinder`
- orthofinder conda package came with the wrong version of diamond, so manually installed diamond and replaced the executable in the conda env directory

~~~
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
cp diamond /hb/home/aanakamo/.conda/envs/orthofinder/bin/
~~~

`conda install -c conda-forge ncbi-datasets-cli`

`conda install -c conda-forge biopython`
