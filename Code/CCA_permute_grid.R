#This function performs CCA.permute over a grid of penalties

CCA_permute_grid <- function(X, Z, chromx = NULL, x_penalty = seq(0,1,0.1), 
z_penalty = seq(0,1,0.1), standardize = FALSE, nperms = 100, 
search_grid = TRUE, filename){

    if(!file.exists(filename)){
        if(search_grid){
            penalty_pairs <- cbind(rep(x_penalty, length(z_penalty)), 
            rep(z_penalty, each = length(x_penalty)))
        }else{
            penalty_pairs <- cbind(x_penalty, z_penalty)
        }
        if(is.null(chromx)){
            perm.results <- CCA.permute(x = X, Z, typex = "standard", typez = "standard", 
            penaltyxs = penalty_pairs[,1], penaltyzs = penalty_pairs[,2], nperms = nperms,
            standardize = standardize)
        }else{
            perm.results <- CCA.permute(x = X, Z, typex = "ordered", typez = "standard", 
            penaltyxs = penalty_pairs[,1], penaltyzs = penalty_pairs[,2], nperms = nperms,
            chromx = chromx, standardize = standardize)
        }
        saveRDS(perm.results, filename)
    }else{
        perm.results <- readRDS(filename)
    }
    return(perm.results)
}