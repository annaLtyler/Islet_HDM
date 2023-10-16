#this function normalizes genoprobs so they can be
#compared to allele loadings. The allele loadings
#are normalized and mean centered as shown in 
#get_allele_loadings()
#ind_id is the IDs of the individuals included
#in the model fitting. This can be derived from
#the expression matrix or the fit models. The
#individuals need to be the same because we 
#mean-center across individuals.

norm_genoprobs <- function(genoprobs, ind_id, verbose = FALSE){

    sub_geno <- lapply(genoprobs, function(x) x[ind_id,,])
    norm_geno <- vector(mode = "list", length = length(sub_geno))
    names(norm_geno) <- names(sub_geno)

    for(chr in 1:length(genoprobs)){
        if(verbose){report.progress(chr, length(genoprobs))}

        n_marker <- dim(sub_geno[[chr]])[3]
        norm_probs = array(NA, dim = dim(sub_geno[[chr]]))
        dimnames(norm_probs) <- dimnames(sub_geno[[chr]])
        
        for(m in 1:n_marker){
            
            # Get current marker
            curr_marker = sub_geno[[chr]][,,m]
            
            # Normalize and center
            norm_factor = 1 / sqrt(rowSums(curr_marker^2))
            norm_marker = scale(diag(norm_factor) %*% curr_marker, center = TRUE, scale = FALSE)
            
            # Multiply by individual loadings to get marker-level loadings
            norm_probs[,,m] = norm_marker
        }

        norm_geno[[chr]] <- norm_probs
    }
    if(verbose){cat("\n")}

    return(norm_geno)

}