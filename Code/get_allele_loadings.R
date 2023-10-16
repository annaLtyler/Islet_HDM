#uses individual-level genome scores and genoprobs
#to calculate allele-level loadings for every marker.

get_allele_loadings <- function(genotype_loadings, genoprobs){

    common.ind <- intersect(rownames(genotype_loadings), rownames(genoprobs[[1]]))

    matched_geno_loadings <- genotype_loadings[common.ind,,drop=FALSE]
    matched_genoprobs <- lapply(genoprobs, function(x) x[common.ind,,])

    allele_loadings = vector(mode = "list", length = length(genoprobs))
    names(allele_loadings) <- names(genoprobs)

    for(chr in 1:length(genoprobs)){
    
        n_marker = dim(genoprobs[[chr]])[3]
        n_allele = dim(genoprobs[[chr]])[2]
    
        chr.loadings = matrix(NA, nrow = n_marker, ncol = n_allele)
        colnames(chr.loadings) = dimnames(genoprobs[[chr]])[[2]]
        rownames(chr.loadings) <- dimnames(genoprobs[[chr]])[[3]]
        
        for(m in 1:n_marker){
            
            # Get current marker
            curr_marker = matched_genoprobs[[chr]][,,m]
            
            # Normalize and center
            norm_factor = 1 / sqrt(rowSums(curr_marker^2))
            norm_marker = scale(diag(norm_factor) %*% curr_marker, center = TRUE, scale = FALSE)
            
            # Multiply by individual loadings to get marker-level loadings
            chr.loadings[m, ] = t(matched_geno_loadings) %*% norm_marker
            
        }
        
        allele_loadings[[chr]] <- chr.loadings
    }

    return(allele_loadings)

}