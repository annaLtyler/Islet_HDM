setwd("~/Documents/Projects/Islets/Islet_HDM")
library(here)
library(pheatmap)

tissue.cols <- as.matrix(read.delim(here("Data", "general", "tissue_colors.txt"), header = FALSE, row.names = 1))[,1]

sep <- readRDS(here("Results", "Paper", "Model_Variance_Explained_tissue_sep-_-complete_mediation.RDS"))
tog <- readRDS(here("Results", "Paper", "Model_Variance_Explained_tissue_together-_-complete_mediation.RDS"))

var.exp.mat <- rbind(unlist(sep), unlist(tog))
rownames(var.exp.mat) <- c("Separate", "Together")

pdf(here("Results", "Paper", "Comparison_Variance_Explained.pdf"), width = 9, height = 5)
layout(matrix(c(1,2), ncol = 2), widths = c(1, 0.3))
par(mar = c(4,4,4,0))
barplot(var.exp.mat, beside = TRUE, density = c(25, 50), col = rep(tissue.cols, each = 2))
mtext("Variance Explained", side = 2, line = 2.5)
par(mar = c(4,0,4,1))
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,1))
legend(x = 0, y = 0.9, fill = "darkgray", density = c(25, 50), legend = c("separate", "together"))

pcor.files <- list.files(here("Results", "Paper"), pattern = "Pcor", full.names = TRUE)
pcor.results <- lapply(pcor.files, readRDS)
pcor.mat <- t(sapply(pcor.results, function(x) x[upper.tri(x)]))
rownames(pcor.mat) <- gsub(".RDS", "", gsub("Pcor_", "", basename(pcor.files)))

barplot(pcor.mat[,c(1,3,2)], beside = TRUE, names = c("Genome-Transcriptome", "Transcriptome-Phenome", "Genome-Phenome"), col = rep(c("darkgray", tissue.cols), 3), ylab = "Correlation", ylim = c(0,0.8))
plot.dim <- par("usr")
dev.off()

pdf(here("Results", "Paper", "Comparison_Correlation.pdf"), width = 9, height = 5)
plot.new()
plot.window(xlim = plot.dim[1:2], ylim = plot.dim[3:4])
abline(h = seq(0, 0.8, 0.1), lty = 2, col = "darkgray")
barplot(pcor.mat[,c(1,3,2)], beside = TRUE, names = c("Genome-Transcriptome", "Transcriptome-Phenome", "Genome-Phenome"), col = rep(c("darkgray", tissue.cols), 3), add = TRUE)
mtext("Correlation", side = 2, line = 2.5)
legend("topright", fill = c("darkgray", tissue.cols), legend = c("All Tissues", names(tissue.cols)), bg = "white")
dev.off()

