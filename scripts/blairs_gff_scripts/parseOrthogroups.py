#!/usr/bin/env python2.7

# This file takes as input an Orthogroups.tsv file (output from OrthoFinder run), a species name, and the column in Orthogroups.tsv 
# pertaining to that species, and outputs a TSV containing a Gene IDs in Column 1 and Orthogroup IDs in Column 2.

# Usage:
# python2 parseOrthogroups.py [path/to/Orthogroups.tsv] [species name of interest from Orthogroups.tsv header] [column number corresponding to species of interest in Orthogroups.tsv]

# Example: 
# python2 parseOrthogroups.py Orthogroups.tsv hSapiens 3

import sys

orthogroups = sys.argv[1]
species_name = sys.argv[2]
column = int(sys.argv[3]) - 1

outfile = orthogroups[:-3] + species_name + '.tsv'

print ""
print "Input file: " + orthogroups
print "Column with info for " + species_name + ": " + sys.argv[3]
print "Saving output to: " + outfile
print ""

gene_to_ortho = {}

with open(orthogroups) as a:
    for line in a.readlines()[1:]:
        line = line.rstrip().split('\t')
        if len(line) > column:
            if len(line[column]) > 0:
                og_id = line[0]
                for entry in line[column].split(","):
                    entry = entry.lstrip().rstrip()
                    if entry not in gene_to_ortho:
                        gene_to_ortho[entry] = og_id


with open(outfile,'w') as out:
    print >> out, '\t'.join(['GeneID','Orthogroup'])    ## changed "Parent" to "GeneID"
    for entry in gene_to_ortho:
        output_entry = [entry.rstrip("_" + species_name),gene_to_ortho[entry]]
        print >> out, '\t'.join(output_entry)

