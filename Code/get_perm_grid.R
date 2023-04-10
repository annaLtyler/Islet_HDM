#This function parses results from CCA_permute_grid

get_perm_grid <- function(cca.perm.results){
    xpen <- cca.perm.results$penaltyxs
    zpen <- cca.perm.results$penaltyzs
    uxpen <- sort(unique(xpen))
    uzpen <- sort(unique(zpen))

    zmat <- cor.mat <- p.mat <- perm.cor.mat <- matrix(NA, nrow = length(uxpen), ncol = length(uzpen))
    rownames(zmat) <- rownames(cor.mat) <- rownames(p.mat) <- rownames(perm.cor.mat) <- uxpen
    colnames(zmat) <- colnames(cor.mat) <- colnames(p.mat) <- colnames(perm.cor.mat) <- uzpen
    for(i in 1:length(xpen)){
        row.idx <- which(uxpen == xpen[i])
        col.idx <- which(uzpen == zpen[i])
        zmat[row.idx, col.idx] <- cca.perm.results$zstats[i]
        cor.mat[row.idx, col.idx] <- cca.perm.results$cors[i]
        p.mat[row.idx, col.idx] <- cca.perm.results$pvals[i]
        perm.cor.mat[row.idx, col.idx] <- mean(cca.perm.results$corperms[i,])
    }
    result <- list("Cor" = cor.mat, "Perm.Cor" = perm.cor.mat, "Z" = zmat, "p" = p.mat)
    return(result)
}
