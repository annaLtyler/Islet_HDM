#This function is based on CCA_permute_grid()
#However, rather than running CCA.permute, 

CCA_permute_components <- function(X, Z, chromx = NULL, x_penalty = seq(0,1,0.1), 
z_penalty = seq(0,1,0.1), nperms = 100, search_grid = TRUE, num_components = NULL,
standardize = FALSE, filename, verbose = TRUE, plot.results = FALSE){

    if(is.null(num_components)){num_components = ncol(Z)}

    if(!file.exists(filename)){
        if(search_grid){
            penalty_pairs <- cbind(rep(x_penalty, length(z_penalty)), 
            rep(z_penalty, each = length(x_penalty)))
        }else{
            penalty_pairs <- cbind(x_penalty, z_penalty)
        }

        observed.cors <- matrix(NA, ncol = num_components, nrow = nrow(penalty_pairs))
        perm.cors <- vector(mode = "list", length = nperms)

        if(is.null(chromx)){
            for(pen in 1:nrow(penalty_pairs)){
                if(verbose){cat("\nRunning CCA for", penalty_pairs[pen,], "\n")}
                observed.result <- CCA(x = X, Z, typex = "standard", typez = "standard", 
                penaltyx = penalty_pairs[pen,1], penaltyz = penalty_pairs[pen,2],
                K = num_components, trace = FALSE)
                observed.cors[pen,] <- observed.result$cors
                
                if(verbose){cat("Running Permutations...\n")}
                perm.cor.mat <- matrix(NA, ncol = num_components, nrow = nperms)
                for(perm in 1:nperms){
                    if(verbose){report.progress(perm, nperms)}
                    rnd.order <- sample(1:nrow(X))
                    perm.result <- CCA(x = X[rnd.order,], Z, typex = "standard", 
                    typez = "standard", penaltyx = penalty_pairs[pen,1], 
                    penaltyz = penalty_pairs[pen,2], K = num_components, 
                    trace = FALSE, standardize = standardize)
                    perm.cor.mat[perm,] <- perm.result$cors
                }
                perm.cors[[pen]] <- perm.cor.mat
            }
        }else{
            #run fused LASSO
            for(pen in 1:nrow(penalty_pairs)){
                if(verbose){cat("Running CCA for", penalty_pairs[pen,], "\n")}
                observed.result <- CCA(x = X, Z, typex = "ordered", typez = "standard", 
                penaltyx = penalty_pairs[pen,1], penaltyz = penalty_pairs[pen,2],
                K = num_components, trace = FALSE, chromx = chromx, standardize = standardize)
                observed.cors[pen,] <- observed.result$cors
                
                if(verbose){"Running Permutations...\n"}
                perm.cor.mat <- matrix(NA, ncol = num_components, nrow = nperms)
                for(perm in 1:nperms){
                    if(verbose){report.progress(perm, nperms)}
                    rnd.order <- sample(1:nrow(X))
                    perm.result <- CCA(x = X[rnd.order,], Z, typex = "ordered", 
                    typez = "standard", penaltyx = penalty_pairs[pen,1], 
                    penaltyz = penalty_pairs[pen,2], K = num_components, trace = FALSE,
                    chromx = chromx, standardize = standardize)
                    perm.cor.mat[perm,] <- perm.result$cors
                }
                perm.cors[[pen]] <- perm.cor.mat
            }
        perm.results <- list("observed.cors" = observed.cors, "perm.cors" = perm.cors)
        saveRDS(perm.results, filename)
        }
    }else{
        perm.results <- readRDS(filename)
    }

    if(plot.results){

        obs.cor <- perm.results$observed.cors
        perm.cor <- perm.results$perm.cors 
        
        layout.mat <- get.layout.mat(nperms)
        layout(layout.mat)
        for(pen in 1:nrow(penalty_pairs)){
            all.cor.mat <- abs(rbind(obs.cor[pen,], perm.cor[[pen]]))
            plot.new()
            plot.window(xlim = c(1, num_components), 
            ylim = c(min(all.cor.mat), max(all.cor.mat)))
            boxplot(abs(perm.cor[[pen]]), add = TRUE)
            points(x = 1:num_components, y = abs(obs.cor[pen,]), col = "red", pch = 16) 
            mtext(side = 3, text = paste(penalty_pairs[pen,], collapse = ", "))
        }
    }

    #calculate p values for each component at each penalty
    obs.cor <- perm.results$observed.cors
    perm.cor <- perm.results$perm.cors 
    pval.mat <- matrix(NA, ncol = num_components, nrow = nrow(penalty_pairs))
    for(pen in 1:nrow(penalty_pairs)){
        pval.mat[pen,] <- sapply(1:ncol(perm.cor[[pen]]), 
        function(x) calc.emp.p(obs.cor[pen,x], perm.cor[[pen]][,x]))
    }

    full.results <- list("observed.cors" = observed.cors, "perm.cors" = perm.cors,
    "penalty_grid" = penalty_pairs, "pvals" = pval.mat)

    return(full.results)
}