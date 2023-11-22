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
- interactive job: `srun -p Instruction -n 8 -N 2 –mem=1G –cpus-per-task=1 –pty /bin/bash`

