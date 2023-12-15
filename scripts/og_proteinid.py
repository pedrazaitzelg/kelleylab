import sys

gff_file = sys.argv[1]
og_tsv = sys.argv[2]

outfile = outfile = og_tsv[:-3] + 'protein_id.tsv'

pID_geneID = {}
with open(gff_file, 'r') as gff:
    for line in gff:
        lst = line.rstrip().split('\t')
        if len(lst) == 9 and lst[2] == "CDS":
            full_att = lst[8].split(';')
            if "protein_id" in line:
                protein_id = [x for x in full_att if "protein_id" in x][0].split('=')[1]
                GeneID = [x for x in full_att if "GeneID" in x][0].split(':')[1].split(',')[0]
                if pID_geneID.get(protein_id):
                    if pID_geneID[protein_id] != GeneID:
                        print("Inconsistent protein_id GeneID pair: " + ", ".join([protein_id, pID_geneID[protein_id], GeneID]))
                else:
                    pID_geneID[protein_id] = GeneID

with open(og_tsv, 'r') as og:
    with open(outfile, 'w') as out:
        for line in og:
            lst = line.rstrip().split()
            if lst[0] == "protein_id":
                print("\t".join(["geneID", "Orthogroup"]))
            else:
                protein_id = lst[0]
                Orthogroup = lst[1]
                geneID = pID_geneID[protein_id]
                out.write("\t".join([geneID, Orthogroup]) + "\n")
