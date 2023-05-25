#this function performs a test for tissue-specific
#eQTL. A pared down version is used in 1a.Tissue_Expression.Rmd
#That version depends on other objects defined in the markdown.
#This version is a standalone version that does a few extra steps 
#for each of use. Because of the extra steps, it takes too long to
#run for use in high volume calculations.
#The arguments are the following:
# transcript.id: the ensembl.id for the transcript of interest.
# tissue.expr: a list containing gene expression for all tissues
# gene.information: a list containing the annot.mrna tables for
# all tissues.
# genoprobs: the genoprobs object for the DO
# K: the kinship object for the DO
# return.data: whether to return only test statistics, or
# additional information about the transcript (see below).
#The function returns a list with the following components:
#if return.data = FALSE
# "tissues.with.expr" = has.expr; # which tissues express the given transcript
# "p.values" = component.p; #p values for effects of tissue, local genotype, and tissue by genotype interaction
#"rel.effects" = relative.effects; # Mean Sq values for all parameters in the model.
#if return.data is TRUE, the following additional elements are returned.
# "genotypes" = nearest.geno; #Genotypes for the marker nearest the gene.
# "expression" = tx.expr; #expression in each of the tissues.
#"chr" = chr.idx; #the chromosome the gene is encoded on
#"gene.name" = gene.name; # the name of the gene matching the transcript ID

get_transcript_effects <- function(transcript.id, tissue.expr, gene.information, 
  genoprobs, K, use.kinship = FALSE, return.data = FALSE){

    tissue.tx <- lapply(tissue.expr, colnames)
    common.ind <- Reduce("intersect", lapply(tissue.expr, rownames))
    ind.matched.expr <- lapply(tissue.expr, function(x) x[common.ind,])
    
    tissue.idx <- lapply(tissue.tx, function(x) which(x == transcript.id))
    has.expr <- which(sapply(tissue.idx, length) > 0)
    
    if(length(has.expr) < 2){return(NULL)}

    tx.expr <- sapply(1:length(has.expr), function(x) ind.matched.expr[[has.expr[[x]]]][,tissue.idx[has.expr][[x]]])
    colnames(tx.expr) <- names(ind.matched.expr)[has.expr]
    
    tx.info  <- gene.tables[[has.expr[1]]][which(gene.tables[[has.expr[1]]]$gene.id == transcript.id),]
    tx.chr <- tx.info$chr
    chr.idx <- which(names(K) == tx.chr)
    gene.name <- tx.info$symbol

    if(!(tx.chr %in% c(1:19, "X"))){return(NULL)}

    nearest.marker <- tx.info$nearest.marker.id
    nearest.geno <- pull_genoprobpos(genoprobs, marker = nearest.marker)

    tx.expr <- lapply(1:length(has.expr), function(x) ind.matched.expr[[has.expr[[x]]]][,tissue.idx[has.expr][[x]]])
    names(tx.expr) <- tissue.names[has.expr]
    long.expr <- unlist(tx.expr)
    tissue.factor <- as.factor(unlist(lapply(1:length(tx.expr), 
      function(x) rep(tissue.names[x], length(tx.expr[[x]])))))
    rep.geno <- Reduce("rbind", lapply(1:length(has.expr), function(x) nearest.geno[rownames(ind.matched.expr[[1]]),]))      
    rep.K <- K[[chr.idx]][rownames(rep.geno), rownames(rep.geno)]
    rownames(rep.K) <- colnames(rep.K) <- names(long.expr)

    if(use.kinship){
      #use cape code to adjust variables
      adj_data <- kin_adjust(kin_obj = rep.K, geno = rep.geno, 
        phenoV = as.matrix(long.expr))
      test <- lm(adj_data$corrected_pheno~tissue.factor*adj_data$corrected_geno)
    }else{
      test <- lm(long.expr~tissue.factor*rep.geno)
      #anova(test)
    }
    test.summary <- anova(test)
    component.p <- test.summary[,"Pr(>F)"][1:3]
    names(component.p) <- rownames(anova(test))[1:3]
    
    
    relative.effects <- test.summary[,"Mean Sq"]
    names(relative.effects) <- c("Tissue", "Genotype", 
      "Interaction", "Residual")

    if(!return.data){
      result = list("tissues.with.expr" = has.expr,
        "p.values" = component.p, "rel.effects" = relative.effects)
    }else{
      result = list("tissues.with.expr" = has.expr,
        "p.values" = component.p, "rel.effects" = relative.effects,
        "genotypes" = nearest.geno, "expression" = tx.expr,
        "chr" = chr.idx, "gene.name" = gene.name)
    }

    return(result)  
}
