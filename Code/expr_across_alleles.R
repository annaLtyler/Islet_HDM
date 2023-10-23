#This is a standalone function that matches that in 1a.Tissue_Expression.Rmd
# It plots expression for all alleles in a given tissue.

expr_across_alleles <- function(transcript.id, tissue.name, tissue.expr, 
    gene.information, genoprobs, K, use.kinship = FALSE){
  
  tx.effects <- get_transcript_effects(transcript.id, tissue.expr, 
    gene.information, genoprobs, K, use.kinship = FALSE, return.data = TRUE)

  if(is.null(tx.effects)){
    plot.text("No comparable expression.")
    return(NULL)
  }

  nearest.geno <- tx.effects$genotypes
  chr.idx <- tx.effects$chr
  expr.list <- tx.effects$expression
  tissue.idx <- which(names(expr.list) == tissue.name)
  relative.effects <- tx.effects$rel.effects

  par(mfrow = c(2,4))
  for(a in 1:ncol(nearest.geno)){
    rounded.geno <- bin.vector(nearest.geno[names(expr.list[[tissue.idx]]),a], c(0, 0.5, 1))
    boxplot(expr.list[[tissue.idx]]~rounded.geno, 
      xlab = "Genotype", ylab = paste(tx.effects$gene.name, "expression"),
      main = names(CCcolors[a]), col = CCcolors[a])
    #plot.with.model(nearest.geno[names(expr.list[[tissue.idx]]),a], 
    #  expr.list[[tissue.idx]], 
    #  xlab = "Genotype", ylab = paste(gene.name, "expression"),
    #  main = names(CCcolors[a]), col = CCcolors[a], report = "cor.test")
  }  
  mtext(tissue.name, side = 3, outer = TRUE, line = -1.5, font = 2)
}
