import sys
import os

os.chdir("/Users/annenakamoto/ROTATION2")   ## locally
sco_file = "Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
og_file = "Orthogroups/Orthogroups.txt"
de_genes_file = "DESEQ2/DE_GENES/ALL.DEgenes.tsv"

### determine number of 1-1 orthologs (SCOs) in the final dataset
SCO = {}
with open(sco_file, 'r') as sco:
    for line in sco:
        og = line.strip()
        SCO[og] = 1

SP_TS_SCO = {}
with open(de_genes_file, 'r') as de:
    for line in de:
        lst = line.strip().split()
        if lst[0] != "species":
            species = lst[0]
            tissue = lst[1]
            og = lst[4]
            sp_ts = (species, tissue)
            if SP_TS_SCO.get(sp_ts):
                SP_TS_SCO[sp_ts][0].append(og)
            else:
                SP_TS_SCO[sp_ts] = [[og], []]
            if SCO.get(og):
                    SP_TS_SCO[sp_ts][1].append(og)

with open("species_tissue_de_scos.tbl.txt", 'w') as out:
    out.write('\t'.join(["species", "tissue", "DE_OG_count", "DE_SCO_count"]) + '\n')
    for sp_ts, ogs in sorted(SP_TS_SCO.items(), key=lambda x: x[0][0]):
        out.write('\t'.join([sp_ts[0], sp_ts[1], str(len(ogs[0])), str(len(ogs[1]))]) + '\n')

with open("species_tissue_de_scos.lst.txt", 'w') as out:
    out.write('\t'.join(["species", "tissue", "DE_OGs", "DE_SCOs"]) + '\n')
    for sp_ts, ogs in sorted(SP_TS_SCO.items(), key=lambda x: x[0][0]):
        out.write('\t'.join([sp_ts[0], sp_ts[1], ";".join(ogs[0]), ";".join(ogs[1])]) + '\n')

### filter out OGs that are only DE in one species (PAV)
### record how many there are
###     Start with the set of all OGs
###     how many are present only in one species and DE in that species (no orthologs in other species)?
###     how many are present in all species and DE in all species?
        
### species, tissue, gene, orthogroup, DE



