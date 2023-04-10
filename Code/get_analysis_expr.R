#This function calculates individual-level CT expression for
#mice in the DO based on a tissue name and CT number.
#you can use imputed expression
#tissue.expr <- readRDS(here("Data", "imputed", "Adjusted_Expression_imputed_local.RDS"))
#tissue.expr <- readRDS(here("Data", "imputed", "Adjusted_Expression_imputed_genetic.RDS"))
#or measured expression
#tissue.expr <- readRDS(here("Data", "Adjusted_Expression.RDS"))
#experiment name refers to which CCA experiment from which you 
#would like to draw the CT weights. Right now there are two options.
#exp.name = "all_traits"
#exp.name = "imputed_local"
#read expression data, either locally imputed or measured
#use either imputed expression
#these files are created by 1.adjust_transcripts_traits.Rmd
#adj.do.expr <- readRDS(here("Data", "imputed", "Adjusted_Expression_DO_imputed_local.RDS"))
#or measured expression
#adj.do.expr <- readRDS(here("Data", "Adjusted_Expression.RDS"))

get_DO_CT_expr <- function(adj.do.expr, exp.name, tissue.name, CT){

    imp.type = imp.type[1]

    #read CT data
    boot.CCA.results <- readRDS(here("Results", "CCA_Clusters", exp.name, "Aggregate.Results.RDS"))
    tissue.idx <- which(names(boot.CCA.results) == tissue.name)
    
    #get weights for the gene expression
    gene.weights <- boot.CCA.results[[tissue.idx]]$u[,CT,drop=FALSE]
    names(gene.weights) <- colnames(adj.do.expr[[tissue.idx]])
    
    all.mouse.id <- unique(unlist(lapply(adj.do.expr, rownames)))
   
    tissue.expr <- adj.do.expr[[which(names(adj.do.expr) == tissue.name)]]
    weighted.genes <- tissue.expr%*%gene.weights
    colnames(weighted.genes) <- paste0(tissue.name, "_CT", CT)
    return(weighted.genes)

}