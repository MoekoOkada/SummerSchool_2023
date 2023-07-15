library(edgeR)
library(DESeq2)
library(dplyr)
library(stringr)
library(cummeRbund)
library(gplot)

# load count data
setwd("/tmp/eagle/work/polyploid-rna-seq-analyses/src/exp")
hal_count <- read.table("hal_counts.tsv", header = T, row.names = 1, sep = "\t")
lyr_count <- read.table("lyr_counts.tsv", header = T, row.names = 1, sep = "\t")

# Make expression matrix with all data
count <- cbind(hal_count, lyr_count)
head(count)

# load gene length of hal
Ahal_gtf <- read.table("../genome/Ahal_v2_2.gtf", header = F, sep = "\t")
Ahal_gtf <- filter(Ahal_gtf, V3 == "transcript")
len <- Ahal_gtf$V5 - Ahal_gtf$V4

# Calculate RPKM
rpkm <- as.data.frame(rpkm(count, len, log = FALSE))
table_rpkm <- as.data.frame(rpkm, n = nrow(count))
head(table_rpkm)
write.table(table_rpkm, file = paste0("homoeolog_RPKM.txt"), col.names = T, row.names = T, sep = "\t")

# make count matrix for each condition (control and zinc)
control <- select(count, -contains("48h"))
head(control)
zinc <- select(count, contains("48h"))
head(zinc)

# set groups
group <- factor(c("kam_hal", "kam_hal", "kam_hal", "kam_lyr", "kam_lyr", "kam_lyr"))
design <- model.matrix(~group)

# DEG detection from zinc treated samples
d1 <- DGEList(counts = zinc, group = group)

# normalize according to UserGuide
d1 <- calcNormFactors(d1)
d1 <- estimateDisp(d1, design)

# plot Multi-dimensional scaling
plotMDS(d1)

# calculate cpn and format it for heatmap
d2 <- cpm(d1, normalized.lib.size = TRUE)
scaledata <- t(scale(t(d2))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata), ]

# Change ids
rownames(scaledata) <- lapply(rownames(scaledata), gsub, pattern = ".t1", replacement = "")

# make a gene set
metal_genes <- read.table("metal_genes_ids.txt", sep = "\t", header = T)
metal_genes[, 4]

# Extract metal genes from deg files
my.geneset <- scaledata[rownames(scaledata) %in% metal_genes[, 4], ]

# clustering
hr <- hclust(as.dist(1 - cor(t(my.geneset), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(my.geneset, method = "spearman")), method = "complete")

# Clustering of metal related genes
png(paste0("Cluster_hm_zinc.png"))
heatmap.2(my.geneset,
  Rowv = as.dendrogram(hr),
  Colv = as.dendrogram(hc),
  col = colorRampPalette(c("blue", "yellow"))(100),
  scale = "row",
  cexCol = 0.7,
  main = "Heatmap.2",
  trace = "none"
)
dev.off()
