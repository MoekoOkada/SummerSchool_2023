library(locfit)
library(doMC)
registerDoMC(12)
source("calcpval.R")

cat("Loading\n") ##
args <- commandArgs(trailingOnly = T)
outfile <- args[1]
rawcount <- read.delim(args[2], header = T)
head(rawcount)
mylabel <- read.table(args[3], header = F, stringsAsFactors = F, sep = ",")
head(mylabel)

MINVALVAR <- 1.0e-5
# MINPSD=1.0e-1
# MINVALVAR=1.0e-5
MINPVAR <- 1.0e-2

cat("Initializing\n") ##
nrow(mylabel)
nRep <- nrow(mylabel)
labelAllRow <- c()
for (i in 1:nRep) {
  for (k in c(1, 2)) {
    mylabel[i, k]
    labelAllRow <- c(labelAllRow, mylabel[i, k])
  }
}
for (i in 1:nRep) {
  for (k in c(3, 4)) {
    labelAllRow <- c(labelAllRow, mylabel[i, k])
  }
}

counts <- rawcount[, labelAllRow]

for (indexCols in list(c(1, 2), c(3, 4))) {
  for (indexRep in 1:nRep) {
    col1 <- mylabel[indexRep, indexCols[1]] # parent1-origin
    col2 <- mylabel[indexRep, indexCols[2]] # parent2-origin
    counts[, col1] <- ifelse(counts[, col1] < counts[, col2] * 0.01, counts[, col2] * 0.01, counts[, col1])
    counts[, col2] <- ifelse(counts[, col2] < counts[, col1] * 0.01, counts[, col1] * 0.01, counts[, col2])
  }
}

OCcolumns <- mylabel[, 1]
ACcolumns <- mylabel[, 2]
ODcolumns <- mylabel[, 3]
ADcolumns <- mylabel[, 4]

ORGcolumns <- c(OCcolumns, ODcolumns)

totalCounts <- counts[, ORGcolumns] + counts[, c(ACcolumns, ADcolumns)]
totalLoggeomeans <- rowMeans(log(totalCounts))
ORGscaledCounts <- exp(log(counts[, ORGcolumns]) - (log(totalCounts) - totalLoggeomeans))
ORGVarList <- apply(ORGscaledCounts, 1, var) * (length(ORGcolumns) - 1) / length(ORGcolumns)
ORGRatioList <- rowMeans(counts[, ORGcolumns] / totalCounts)
ORGRatioVarList <- apply(counts[, ORGcolumns] / totalCounts, 1, var) * (length(ORGcolumns) - 1) / length(ORGcolumns)
ORGRatioSDList <- sqrt(ORGRatioVarList)

ctrlCounts <- counts[, OCcolumns] + counts[, ACcolumns]
ctrlLoggeomeans <- rowMeans(log(ctrlCounts))
OCscaledCounts <- exp(log(counts[, OCcolumns]) - (log(ctrlCounts) - ctrlLoggeomeans))
OCVarList <- apply(OCscaledCounts, 1, var) * (length(OCcolumns) - 1) / length(OCcolumns)
OCRatioList <- rowMeans(counts[, OCcolumns] / ctrlCounts)
OCRatioVarList <- apply(counts[, OCcolumns] / ctrlCounts, 1, var) * (length(OCcolumns) - 1) / length(OCcolumns)

objCounts <- counts[, ODcolumns] + counts[, ADcolumns]
objLoggeomeans <- rowMeans(log(objCounts))
ODscaledCounts <- exp(log(counts[, ODcolumns]) - (log(objCounts) - objLoggeomeans))
ODVarList <- apply(ODscaledCounts, 1, var) * (length(ODcolumns) - 1) / length(ODcolumns)
ODRatioList <- rowMeans(counts[, ODcolumns] / objCounts)
ODRatioVarList <- apply(counts[, ODcolumns] / objCounts, 1, var) * (length(ODcolumns) - 1) / length(ODcolumns)

cat("Fitting\n")
finiteval <- is.finite(ORGVarList) & totalLoggeomeans > 0 & ORGVarList > 1.0e-10
fit <- locfit(log(ORGVarList[finiteval]) ~ lp(totalLoggeomeans[finiteval], ORGRatioList[finiteval], scale = T), family = "gaussian", maxk = 200)

geneList <- c()
step <- 200
lastnum <- nrow(rawcount)
sampleList <- 1:lastnum
repeatnum <- lastnum %/% step
rest <- lastnum %% step
if (rest > 0) repeatnum <- repeatnum + 1

# dopar
pvalListList <- foreach(i = 1:repeatnum) %dopar% { # length(sampleList)
  firstnum <- (i - 1) * step + 1
  termnum <- min(i * step, lastnum)
  cat(
    "Processing i=", i, "/", repeatnum,
    " from ", firstnum, " to ", termnum, "\n"
  )
  sigChangeMH(firstnum:termnum, numSample = 10000, verbose = FALSE)
}
pvalList <- unlist(pvalListList)
gene <- as.character(rawcount[sampleList, "gene"])
ctrlFstCounts <- rowSums(rawcount[sampleList, OCcolumns])
ctrlSndCounts <- rowSums(rawcount[sampleList, ACcolumns])
objFstCounts <- rowSums(rawcount[sampleList, ODcolumns])
objSndCounts <- rowSums(rawcount[sampleList, ADcolumns])

cat("Perform Fisher's exact test\n")
pfisher <- sapply(1:length(sampleList), function(i) {
  fisher.test(
    matrix(
      c(ctrlFstCounts[i], ctrlSndCounts[i], objFstCounts[i], objSndCounts[i]),
      nrow = 2
    )
  )$p
})

cat("Perform Chisqr Test\n")

pchisqr <- sapply(1:length(sampleList), function(i) {
  ifelse(
    ctrlFstCounts[i] == 0 & ctrlSndCounts[i] == 0 &
      objFstCounts[i] == 0 & objSndCounts[i] == 0,
    1.0,
    chisq.test(
      matrix(
        c(
          ctrlFstCounts[i],
          ctrlSndCounts[i],
          objFstCounts[i],
          objSndCounts[i]
        ),
        nrow = 2
      )
    )$p.value
  )
})

# pvalRes["chisqr"] <- pchisqr
# pvalRes["chisqr.adj"] <- p.adjust(pchisqr)

pvalRes <- data.frame(
  gene = gene, pval = pvalList, padj = p.adjust(pvalList, "BH"),
  ctrlFirst = ctrlFstCounts, ctrlSecond = ctrlSndCounts,
  objFirst = objFstCounts, objSecond = objSndCounts,
  ctrlRatio = ctrlFstCounts / (ctrlFstCounts + ctrlSndCounts),
  objRatio = objFstCounts / (objFstCounts + objSndCounts),
  ratioSD = ORGRatioSDList[sampleList],
  fisher = pfisher, fisher.adj = p.adjust(pfisher, "BH"),
  chisqr = pchisqr, chisqr.adj = p.adjust(pchisqr, "BH")
)

write.table(pvalRes, file = outfile, sep = "\t", quote = F, row.names = F)
