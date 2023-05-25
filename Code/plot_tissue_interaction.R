#This function is a standalone version of a function of the
#same name in 1a.Tissue_Expression.Rmd. It requires more inputs
#and does more calculations than that function, but can operate
#on its own. It takes the same inputs as get_transcript_effects()
#plus text.gap and text.shift for use in barplot_with_num()

plot_tissue_interaction <- function(transcript.id, tissue.expr, 
    gene.information, genoprobs, K, use.kinship = FALSE, 
    text.gap = 0.1, text.shift = 0.05){

    tx.effects <- get_transcript_effects(transcript.id = transcript.id, 
        tissue.expr = tissue.expr, gene.information = gene.information, 
        genoprobs = genoprobs, K = K, use.kinship = use.kinship, 
        return.data = TRUE)

    if(is.null(tx.effects)){
      return("Not expressed in multiple tissues")
    }

    has.expr <- tx.effects$tissues.with.expr
    nearest.geno <- tx.effects$genotypes
    chr.idx <- tx.effects$chr
    expr.mat <- Reduce("cbind", tx.effects$expression)
    colnames(expr.mat) <- names(has.expr)
    relative.effects <- tx.effects$rel.effects
    total.effect <- sum(relative.effects)
    prop.effect <- relative.effects/total.effect

    test <- apply(expr.mat, 2, function(x) fit1(nearest.geno, x, kinship = K[[chr.idx]]))
    tissue.lod <- sapply(test, function(x) x$lod)

    layout(matrix(c(1,2,3,4,4,4), nrow = 2, byrow = TRUE))
    par(mar = c(4,4,4,4))
    boxplot(expr.mat, main = "Expression by tissue", col = tissue.cols[has.expr], 
      las = 2, ylab = "TPM")
    bar.y <- get_plot_bounds(0, round(max(tissue.lod)), scale.factor = 0.2)
    barplot_with_num(round(tissue.lod), main = "LOD score by tissue", 
      col = tissue.cols[has.expr], las = 2, text.gap = text.gap, 
      text.shift = text.shift)
    barplot_with_num(signif(prop.effect, 2), main = "Relative Effects",
      ylab = "Mean Effect")
    
    all.coef <- sapply(test, function(x) x$coef)
    #pheatmap(all.coef[1:8,], scale = "none")
    par(mar = c(4,4,0,4))
    a <- barplot(t(all.coef[1:8,]), beside = TRUE, col = tissue.cols[has.expr], 
      ylab = "Allele Effects", names = rep("", 8))
    bar.range <- max(all.coef[1:8,]) - min(all.coef[1:8,])
    par(xpd = NA)
    text(x = colMeans(a), y = rep(min(all.coef[1:8,])-bar.range*0.05, 8), 
      labels = names(CCcolors), col = CCcolors, font = 2, cex = 2)
    par(xpd = FALSE)

    #I think this barplot is less useful
    #barplot(all.coef[1:8,], beside = TRUE, col = CCcolors, ylab = "Allele Effects")
    
    gene.name <- tx.effects$gene.name
    int.p <- signif(tx.effects$p.values[3], 2)
    mtext(paste0(gene.name, " (p = ", int.p, ")"), 
      side = 3, outer = TRUE, line = -1.5, font = 4, cex = 1.5)
 
}
