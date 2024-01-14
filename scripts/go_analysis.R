install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

library("optparse")
library(stringr)
option_list = list(
  make_option(c("-s", "--subset"), type="character", default=NULL, 
              help="species and tissue [REQUIRED]", metavar="character"),
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
subset = opt$subset

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

setwd("/hb/scratch/aanakamo/kelleylab_rotation/go_analysis")

# reading in data from deseq2
df = read.table("ALL.DEgenes.with_human.tsv", header=TRUE)
colnames(df)
df <- subset(df, df$human_RefSeq_id != "None")   ### remove rows with no human_RefSeq_id
df

### subsets
adipose.brown_bear <- subset(df, df$species == "brown_bear" & df$tissue == "adipose")
dim(adipose.brown_bear)
adipose.syrian_hamster <- subset(df, df$species == "syrian_hamster" & df$tissue == "white_adipose")
dim(adipose.syrian_hamster)
adipose.dwarf_lemur <- subset(df, df$species == "dwarf_lemur" & df$tissue == "white_adipose")
dim(adipose.dwarf_lemur)

liver.brown_bear <- subset(df, df$species == "brown_bear" & df$tissue == "liver")
dim(liver.brown_bear)
liver.brandts_bat <- subset(df, df$species == "brandts_bat" & df$tissue == "liver")
dim(liver.brandts_bat)

muscle.brown_bear <- subset(df, df$species == "brown_bear" & df$tissue == "muscle")
dim(muscle.brown_bear)
muscle.arctic_ground_squirrel <- subset(df, df$species == "arctic_ground_squirrel" & df$tissue == "muscle")
dim(muscle.arctic_ground_squirrel)

brain.13_lined_ground_squirrel <- subset(df, df$species == "13_lined_ground_squirrel")
dim(brain.13_lined_ground_squirrel)
brain.brandts_bat <- subset(df, df$species == "brandts_bat" & df$tissue == "brain")
dim(brain.brandts_bat)

# setup gene_list
tissue.species <- eval(as.name(subset))
original_gene_list <- tissue.species$log2FoldChange   # we want the log2 fold change 
names(original_gene_list) <- tissue.species$human_RefSeq_id   # name the vector
gene_list<-na.omit(original_gene_list)   # omit any NA values 
gene_list <- sort(gene_list, decreasing = TRUE)   # sort the list in decreasing order (required for clusterProfiler)

# gsea
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ACCNUM", 
             nPerm = 100, 
             minGSSize = 10, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

write.csv(gse, file = paste0(subset, ".gse_result.csv") row.names = FALSE)

# dotplot
require(DOSE)
png(paste0(subset, ".dotplot.png"), width = 6, height = 6, units = "in", res = 300)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off() 
