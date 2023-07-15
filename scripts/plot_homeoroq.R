MINRPKM <- 0.2 # 10^-1
MINRATIOSD <- 0.3
PADJ <- 0.05

args <- commandArgs(trailingOnly = T)
# rpkm_csv <- args[1]
pval_mean_xls <- args[1]
control <- args[2]
target <- args[3]
#  e.g.
#  Rscript plot.R rpkms.csv pval_mean.xls w2221_L0 w2221_L48
# rpkm_raw <- read.table(rpkm_csv, header=T, sep=",")
sigchange <- read.table(pval_mean_xls, header = T, sep = "\t")
# rpkm <- merge(rpkm_raw, sigchange[,c("gene", "pval", "padj")], by="gene")

#  IRF0418_vs_IRF0516
# control <- "IRF0418"
# target <- "IRF0516"
# ctrlHalRPKM <- rpkm[,"RPKM_hal.ctrl_all"]
# ctrlLyrRPKM <- rpkm[,"RPKM_lyr.ctrl_all"]
# coldHalRPKM <- rpkm[,"RPKM_hal.cold_all"]
# coldLyrRPKM <- rpkm[,"RPKM_lyr.cold_all"]


# FRbrF0418_1_cama_genome_orig
#  FRbrF0418_1_cama_genome_other
#  FRslF0418_1_cama_genome_orig
#  FRslF0418_1_cama_genome_other

# ctrlHalRPKM <- rpkm[,"FRbrF0418_orig_all"]
# ctrlLyrRPKM <- rpkm[,"FRbrF0418_other_all"]
# coldHalRPKM <- rpkm[,"FRslF0418_orig_all"]
# coldLyrRPKM <- rpkm[,"FRslF0418_other_all"]

# ctrlHalRPKM <- rpkm[, paste(control, "_orig_all", sep="")]
# ctrlLyrRPKM <- rpkm[, paste(control, "_other_all", sep="")]
# coldHalRPKM <- rpkm[, paste(target, "_orig_all", sep="")]
# coldLyrRPKM <- rpkm[, paste(target, "_other_all", sep="")]

ctrlHalRPKM <- sigchange[, "ctrlFirst"]
ctrlLyrRPKM <- sigchange[, "ctrlSecond"]
coldHalRPKM <- sigchange[, "objFirst"]
coldLyrRPKM <- sigchange[, "objSecond"]


ctrlRatio <- ctrlHalRPKM / (ctrlHalRPKM + ctrlLyrRPKM)
hasTagInCtrl <- (ctrlHalRPKM + ctrlLyrRPKM) >= MINRPKM
coldRatio <- coldHalRPKM / (coldHalRPKM + coldLyrRPKM)
hasTagInCold <- (coldHalRPKM + coldLyrRPKM) >= MINRPKM

# sigGenes <- rep(FALSE, nrow(rpkm))
sigGenes <- rep(FALSE, nrow(sigchange))
if (TRUE) {
     # sigGenes <- rpkm[,"padj"] < 0.05
     sigGenes <- sigchange[, "padj"] < PADJ
}
ratiosdFiltering <- sigchange[, "ratioSD"] < MINRATIOSD

pdf(paste("ratio_overdispersion_", control, "_vs_", target, ".pdf", sep = ""))
r_col <- cor(ctrlRatio[hasTagInCtrl & hasTagInCold], coldRatio[hasTagInCtrl & hasTagInCold], method = "pearson")
plot(ctrlRatio[hasTagInCtrl & hasTagInCold],
     coldRatio[hasTagInCtrl & hasTagInCold],
     col = "#55555588",
     pch = 20,
     xlim = c(0, 1), ylim = c(0, 1),
     # xlab="A-origin ratio in FRbrF0418", ylab="A-origin ratio in FRslF0418", cex=1.0) # change "cold" to "stress"
     main = paste(sum(hasTagInCtrl & hasTagInCold), "genes, Sig:", sum(hasTagInCtrl & hasTagInCold & sigGenes & ratiosdFiltering), "genes (padj <", PADJ, "),", sprintf("R=%.2f", r_col), ",\n TPM >", MINRPKM, ", ratioSD <", MINRATIOSD),
     xlab = paste("A-origin ratio in", control), ylab = paste("A-origin ratio in", target), cex = 1.0
) # change "cold" to "stress"
par(new = T)
plot(ctrlRatio[hasTagInCtrl & hasTagInCold & sigGenes & ratiosdFiltering],
     coldRatio[hasTagInCtrl & hasTagInCold & sigGenes & ratiosdFiltering],
     col = "#AA000077",
     pch = 20,
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "", ylab = "", cex = 1.5
)
cat("Points in all:", sum(hasTagInCtrl & hasTagInCold), "genes\n")
cat("Signif:", sum(hasTagInCtrl & hasTagInCold & sigGenes & ratiosdFiltering), "genes\n")
dev.off()

# r_col <- cor(ctrlRatio[hasTagInCtrl&hasTagInCold], coldRatio[hasTagInCtrl&hasTagInCold], method="pearson")
cat("R=", r_col, "\n")

pdf(paste("hist_overdispersion_", control, "_vs_", target, ".pdf", sep = ""))
par(mfrow = c(2, 1))
hist(ctrlRatio[hasTagInCtrl & hasTagInCold], xlab = paste("A-origin ratio in", control), main = "", col = "gray", yaxt = "n")
axis(side = 2, at = seq(0, 2500, 500), labels = c(0, 500, "", 1500, "", 2500), las = 2)
hist(coldRatio[hasTagInCtrl & hasTagInCold], xlab = paste("A-origin ratio in", target), main = "", col = "gray", yaxt = "n")
axis(side = 2, at = seq(0, 2500, 500), labels = c(0, 500, "", 1500, "", 2500), las = 2)
dev.off()
