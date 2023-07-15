# Command lists 2

Date: 17-Jul-2023  
Author: Moeko Okada  

List of commands used in the polyploidy summer school.  

## Part 2

6. Expression analysis
   1. Clustering
   2. Homoeolog ratio test
   3. GO enrichment analysis (if we have time...)

## 5. Download the whole expression data

```bash
$ cd src
src $ wget -O exp.tar.gz https://drive.switch.ch/index.php/s/4JaGxRe5nJvJ8WX/download
src $ ls
1_qc  2_map  3_eagle  4_count  exp.tar.gz
fastq  fastq.tar.gz  genome  genome.tar.gz  scripts
src $ tar xvzf exp.tar.gz
src $ ls
1_qc  2_map  3_eagle  4_count  exp  exp.tar.gz
fastq  fastq.tar.gz  genome  genome.tar.gz  scripts
src $ ls exp
hal_counts.tsv      lyr_counts.tsv      metal_genes_ids.txt
src $
```

## 6. Expression analysis

Use R for calculation nad visualization.

Please move to the R console.

<br>

Required package information is written in [install.R](https://github.com/MoekoOkada/SummerSchool_2023/blob/main/scripts/install.r).


### Load Required libraries

```r
library(edgeR) # normalisation, deg detection
library(DESeq2) # normalisation, deg detection
library(dplyr) # data formatting
library(stringr) # data formatting
library(gplot) # plot heatmap
```

### Load data

```r
## load count data
# set working directory
setwd("/tmp/eagle/work/polyploid-rna-seq-analyses/src/exp")

# load count data of halleri side
hal_count <- read.table("hal_counts.tsv", header = T, row.names = 1, sep = "\t")

# check hal count data
head(hal_count)

# load count data of lyrata side
lyr_count <- read.table("lyr_counts.tsv", header = T, row.names = 1, sep = "\t")

# check lyr count data
head(lyr_counr)
```

### Data prep. for EAGLE-RC

1. Make whole count matrix
2. Load the gene length info of halleri
3. calculate RPKM using DESeq2

```R
# bind count data
count <- cbind(hal_count, lyr_count) # cbind merges the two data frame based on rownames
# check data
head(count)

# load gene length of hal
Ahal_gtf <- read.table("../genome/Ahal_v2_2.gtf", header = F, sep = "\t")
Ahal_gtf <- filter(Ahal_gtf, V3 == "transcript") # extract transcript information
len <- Ahal_gtf$V5 - Ahal_gtf$V4 + 1 # define length based on start and end position of each genes

# Calculate RPKM
rpkm <- as.data.frame(rpkm(count, len, log = FALSE)) # rpkm function in DESeq2 package
table_rpkm <- as.data.frame(rpkm, n = nrow(count)) # convert to table
head(table_rpkm)
write.table(table_rpkm, file = paste0("homoeolog_RPKM.txt"), col.names = T, row.names = T, sep = "\t") # Output
```

### make count matrix for each condition (control and zinc) for expression analysis

```R
# select function in dplyr package
control <- select(count, -contains("48h")) # Select columns without "48h"
head(control)
zinc <- select(count, contains("48h")) # Select columns with "48h"
head(zinc)
```

### Expression analysis: two group comparison

Use edgeR to compare halleri side and kamchatica side under zinc treatment.

```R
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
```

### Heatmap of metal-related gene

First, format data to draw heatmap.

```R
# calculate cpm and format it for heatmap
d2 <- cpm(d1, normalized.lib.size = TRUE)
scaledata <- t(scale(t(d2))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata), ]

# Change ids
rownames(scaledata) <- lapply(rownames(scaledata), gsub, pattern = ".t1", replacement = "")
```

### Heatmap of metal-related gene

Extract metal-related gene from expression table.

```R
# make a gene set
metal_genes <- read.table("metal_genes_ids.txt", sep = "\t", header = T)
metal_genes[, 4]

# Extract metal genes from deg files
my.geneset <- scaledata[rownames(scaledata) %in% metal_genes[, 4], ]
```

### Heatmap of metal-related gene

Hierarchical clustering

```R
# clustering
hr <- hclust(as.dist(1 - cor(t(my.geneset), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(my.geneset, method = "spearman")), method = "complete")
```

### Heatmap of metal-related gene

Plot heatmap based on the clustering information.

Use `heatmap.2` function in `gplot` package.

```R
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
```

## Homoeologous ratio test

### Make sure you are in `scripts` directory.

```bash
$ cd src/scripts
scripts $ R --vanilla --slave --args ../exp/pval.txt ../exp/homoeolog_RPKM.txt label.txt < calcpval_one.R
```