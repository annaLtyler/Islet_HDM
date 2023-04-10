get_standard_gene_info <- function(values, filter = "ensembl_gene_id"){

    require(biomaRt)
    mus <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") 

    gene.info <- getBM(filters = filter, attributes = c("ensembl_gene_id", 
    "external_gene_name", "chromosome_name", "start_position", "end_position"), 
    values = values, mart = mus)

    return(gene.info)
}