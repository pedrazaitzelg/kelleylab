import sys
import os

os.chdir("/Users/annenakamoto/ROTATION2/DESEQ2/GO_ANALYSIS")    ## locally
map_file = "Orthogroups_ALL.protein_id.tsv"     ## Cols: "species, protein_id, GeneID, Orthogroup"
og_file = "Orthogroups_withHuman.txt"           ## Format: "OG#: Prot1 Prot2 etc"
de_file = "ALL.DEgenes.tsv"                     ## Cols: "species, tissue, gene_name, gene_id, Orthogroup, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, IHW_pvalue"
out_file = "ALL.DEgenes.with_human.tsv"         ## Cols: "species, tissue, gene_name, gene_id, Orthogroup, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, IHW_pvalue, human_RefSeq_id"

### using map_file, store mapping of (species, RefSeq_id) to (species, gene_id)
TO_REFSEQ = {}
with open(map_file, 'r') as m:
    for line in m:
        if "species" not in line:
            lst = line.strip().split()
            species = lst[0]; gene_id = lst[2]; RefSeq_id = lst[1]
            TO_REFSEQ[(species, RefSeq_id)] = (species, gene_id)

### using og_file, store mapping of (species, gene_id) to human_RefSeq_id
### if more than one human protein in an OG, just choose one (arbitrarily choose the first one)
TO_HUMAN_REFSEQ = {}
with open(og_file, 'r') as f:
    for line in f:
        s = line.strip().split(":")
        og = s[0]
        prot_lst = s[1].strip().split()
        human_prots = [p for p in prot_lst if len(p.split("_")) <= 2]
        other_prots = [p for p in prot_lst if len(p.split("_")) > 2]
        if len(human_prots) > 0:
            rep_human = human_prots[0]
        else:
            rep_human = "None"
        for prot in other_prots:
            refseq = "_".join(prot.split("_")[:2])
            species = "_".join(prot.split("_")[2:])
            TO_HUMAN_REFSEQ[TO_REFSEQ[(species, refseq)]] = rep_human

### add human_RefSeq_id column to the de_file and write to out_file
with open(de_file, 'r') as de:
    with open(out_file, 'w') as out:
        for line in de:
            lst = line.strip().split()
            if lst[0] == "species":
                out.write("\t".join(lst + ["human_RefSeq_id"]) + "\n")
            else:
                species = lst[0]
                gene_id = lst[3]
                human_refseq = "None"
                if TO_HUMAN_REFSEQ.get((species, gene_id)):
                    human_refseq = TO_HUMAN_REFSEQ[(species, gene_id)]
                out.write("\t".join(lst + [human_refseq]) + "\n")
