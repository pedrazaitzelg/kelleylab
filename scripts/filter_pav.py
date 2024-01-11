import sys
import os

os.chdir("/Users/annenakamoto/ROTATION2")   ## locally
sco_file = "Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
og_file = "Orthogroups/Orthogroups.txt"
de_genes_file = "DESEQ2/DE_GENES/ALL.DEgenes.tsv"
ogs_present_in_all_sp = "og_size_dist_10.txt"

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

### separate OGs into two sets:
###     present in all species
###     absent from some species

PRES_ALL = {}        
with open(ogs_present_in_all_sp, 'r') as f:
    for line in f:
        lst = line.strip().split()
        PRES_ALL[lst[0]] = 1

with open(de_genes_file, 'r') as de:
    with open("DESEQ2/DE_GENES/ALL.DEgenes.noZERO.tsv", 'w') as noZ:
        with open("DESEQ2/DE_GENES/ALL.DEgenes.ZEROs.tsv", 'w') as Zs:
            for line in de:
                lst = line.strip().split()
                og = lst[4]
                if og == "Orthogroup":
                    noZ.write(line)
                    Zs.write(line)
                elif PRES_ALL.get(og):
                    noZ.write(line)
                else:
                    Zs.write(line)
