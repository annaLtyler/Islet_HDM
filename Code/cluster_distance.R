#This function takes in a square matrix (probably a correlation)
#matrix and a vector of cluster IDs. It calculates the average 
#value within each cluster, and between cluster pairs.

cluster_distance <- function(mat, cID, plot.results = FALSE, ylim = NULL){

    u_clust <- sort(unique(cID))
    n_clust  <- length(u_clust)
    dist.mat <- matrix(NA, n_clust, n_clust)
    rownames(dist.mat) <- colnames(dist.mat) <- u_clust

    clust_pairs <- pair.matrix(u_clust, self.pairs = TRUE)
    all.vals <- vector(mode = "list", length = nrow(clust_pairs))
    names(all.vals) <- apply(clust_pairs, 1, function(x) paste(x, collapse = "_"))

    for(i in 1:nrow(clust_pairs)){
        c1 <- clust_pairs[i,1]
        c2 <- clust_pairs[i,2]
        c1.idx <- which(cID == c1)
        c2.idx <- which(cID == c2)

        #c.order <- order(cID)
        #test <- mat
        #test[c1.idx, c2.idx] <- NA
        #pheatmap(test[c.order, c.order], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)
        c.vals <- mat[c1.idx, c2.idx]
        all.vals[[i]] <- c.vals
        c.dist <- mean(c.vals, na.rm = TRUE)
        dist.mat[c1,c2] <- dist.mat[c2,c1] <- c.dist
    }

    #imageWithText(round(dist.mat, 2), split.at.vals = TRUE, col.scale = c("blue", "brown"), grad.dir = "ends")

    if(plot.results){
        #within.vals <- diag(dist.mat)
        #between.vals <- dist.mat[upper.tri(dist.mat, diag = FALSE)]
        #boxplot(list("Within-Cluster" = within.vals, "Between-Cluster" = between.vals),
        #    ylim = ylim)
        #abline(h = 0)

        boxplot(all.vals)
        on.diag <- which(apply(clust_pairs, 1, function(x) x[1] == x[2]))
        off.diag <- which(apply(clust_pairs, 1, function(x) x[1] != x[2]))
        boxplot(list("Within-Cluster" = unlist(all.vals[on.diag]), 
            "Between-Cluster" = unlist(all.vals[off.diag])))
        abline(h = 0)

    }

     return(dist.mat)

}
