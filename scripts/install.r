# This file contains packages which should be added to the notebook
# during the build process. It is standard R code which is run during
# the build process and typically comprises a set of `install.packages()`
# commands.
#
# For example, remove the comment from the line below if you wish to
# install the `ggplot2` package.
#

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cranPackages <- c(
  "tidyverse",
  "plotly",
  "pheatmap",
  "patchwork",
  "Seurat",
  "dplyr",
  "locfit",
  "doMC"
)
install.packages(cranPackages)

biocPackages <- c(
  "DESeq2",
  "cummeRbund",
  "edgeR",
  "tximport",
  "topGO",
  "ggplot2"
)
BiocManager::install(biocPackages)
