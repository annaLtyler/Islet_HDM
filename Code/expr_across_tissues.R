
#This is a standalone function similar to one in 1a.Tissue_Expression.Rmd
# It plots expression across tissues for a given allele.


expr_across_tissues <- function(transcript.id, allele.name, tissue.expr, 
    gene.information, genoprobs, K, use.kinship = FALSE){

  allele.idx <- c(which(names(CCcolors) == allele.name), which(LETTERS[1:8] == allele.name))
  allele.label <- names(CCcolors)[allele.idx]

  tx.effects <- get_transcript_effects(transcript.id, tissue.expr, 
    gene.information, genoprobs, K, use.kinship = FALSE, return.data = TRUE)
  
  if(is.null(tx.effects)){
    plot.text("No comparable expression.")
    return(NULL)
  }

  has.expr <- tx.effects$tissues.with.expr
  nearest.geno <- tx.effects$genotypes
  chr.idx <- tx.effects$chr
  expr.list <- tx.effects$expression
  relative.effects <- tx.effects$rel.effects  
  gene.name <- tx.effects$gene.name

  layout(get.layout.mat(length(expr.list)))
  for(tx in 1:length(expr.list)){
    rounded.geno <- bin.vector(nearest.geno[names(expr.list[[tx]]),allele.idx], c(0, 0.5, 1))
    boxplot(expr.list[[tx]]~rounded.geno, col = tissue.cols[has.expr[tx]],
      xlab = "Genotype", ylab = paste(gene.name, "expression"),
      main = names(expr.list)[tx])
      #plot.with.model(nearest.geno[names(expr.list[[tx]]),allele.idx], expr.list[[tx]], 
      #xlab = "Genotype", ylab = paste(gene.name, "expression"),
      #main = names(expr.list)[tx], col = tissue.cols[has.expr[tx]],
      #report = "cor.test")
  }  
  mtext(allele.label, side = 3, outer = TRUE, line = -1.5, font = 2)

}