### 2014/08/11, fix mereging

pvalmean <- function(infilenames, outfilename) {
  cat(1, "\n")
  curtable <- read.table(infilenames[1], header = T)
  ptable <- curtable[, 1:2]
  for (i in 2:length(infilenames)) {
    cat(i, "\n")
    curtable <- read.table(infilenames[i], header = T)
    ptable <- merge(ptable, curtable[, 1:2], by = "gene")
  }

  ptable$pval <- rowMeans(ptable[2:ncol(ptable)])
  ptable$padj <- p.adjust(ptable$pval, "BH")
  ptable <- ptable[order(ptable$gene), ]
  curtable <- curtable[order(curtable$gene), ]
  curtable$pval <- ptable$pval
  curtable$padj <- ptable$padj

  write.table(curtable, file = outfilename, row.names = F, sep = "\t", quote = F)
}
