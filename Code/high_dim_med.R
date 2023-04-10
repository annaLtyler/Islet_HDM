#This function performs high-dimensional mediation
#with either kernelized matrices, or data matrices.
#all matrices must have rownames that can be aligned
#by the function
#min.weight.diff is the 
#kernel.c, kernel.m, and kernel.o indicate whether the causal,
#mediating, and outcome matrices are kernelized


high_dim_med <- function(causal.matrix, mediating.matrix, outcome.matrix, 
    min.weight.diff = 1e-3, max.iter = 15, 
    scheme = c("centroid", "horst", "factorial"), verbose = FALSE, 
    kernel.c = TRUE, kernel.m = TRUE, kernel.o = TRUE, use.partial.cor = TRUE){

    scheme <- scheme[1]

    common.ind <- Reduce("intersect", list(rownames(causal.matrix), 
        rownames(mediating.matrix), rownames(outcome.matrix)))

    if(kernel.c){
        g <- causal.matrix[common.ind, common.ind]
    }else{
        g <- causal.matrix[common.ind,,drop=FALSE]
    }

    if(kernel.m){
        t <- mediating.matrix[common.ind, common.ind]
    }else{
        t <- mediating.matrix[common.ind,,drop=FALSE]
    }

    if(kernel.o){
        p <- outcome.matrix[common.ind, common.ind]
    }else{
        p <- outcome.matrix[common.ind,,drop=FALSE]
    }

    check_stop <- function(initial_weights, curr_weights, last_diff){
       
        weight.diff <- abs(initial_weights - curr_weights)

        #check to see if we are converging. 
        #The new weight diffs should be 
        #smaller than the last ones.
        converging <- all(weight.diff < last_diff)

        #check to see if we have reached our minimum 
        #change in weight to meet our stopping criterion.
        reached.min <- all(weight.diff <= min.weight.diff)

        result <- list("reached.min" = reached.min, 
                    "converging" = converging,
                    "weight_diff" = weight.diff)
        return(result)
    }

    decide_stop <- function(stopping_criteria, n.iter){
        do.we.stop = FALSE
        reason <- "keep going"

        #if we have met the minimum weight change criterion, then stop
        if(stopping_criteria[[1]]){
            do.we.stop = TRUE
            reason <- "Reached minimum weight change."
            if(verbose){cat(reason, "\n")}
            
        }

        #if we are not converging, then stop
        if(!stopping_criteria[[2]]){
            do.we.stop = TRUE
            reason <- "Not converging"
            if(verbose){cat(reason, "\n")}
        }

        if(n.iter >= max.iter){
            do.we.stop = TRUE
            reason <- "Reached maximum number of iterations."
            if(verbose){cat(reason, "\n")}
        }

    result <- list("do.we.stop" = do.we.stop, "reason" = reason)
    return(result)
    }

    check_signs <- function(curr_scores){
        if(use.partial.cor){
            score_cor <- pcor.shrink(curr_scores, verbose = FALSE)
        }else{
            score_cor <- cor(curr_scores)
        }
        #cor(curr_scores)
        xm <- score_cor[1,2]
        my <- score_cor[2,3]
        xy <- score_cor[1,3]
        flag <- "coherent"

        if(sign(xm) == -1 && sign(my) == 1){
            curr_scores[,1] <- curr_scores[,1] * -1
            if(sign(xy) == 1){flag = "frustrated"}
        }

        if(sign(xm) == -1 && sign(my) == -1){
            curr_scores[,2] <- curr_scores[,2] * -1
            if(sign(xy) == -1){flag = "frustrated"}
        }

        if(sign(xm) == 1 && sign(my) == -1){
            curr_scores[,3] <- curr_scores[,3] * -1
            if(sign(xy) == 1){flag == "frustrated"}
        }
        
        result <- list("scores" = curr_scores, "flag" = flag)
        return(result)
    }

    A = list(g = g, t = t, p = p)

    weight.mat = 0.5 * matrix(c(0,0,0,1,0,0,0,1,0), 3, 3)    
    weight.mat = weight.mat + t(weight.mat)

    
    #set initial conditions for stopping criteria checks
    stop.now <- FALSE
    W1 <- W2 <- 0.5
    last_diff <- c(Inf, Inf)
    iter = 1
    # Loop over "EM" iterations
    while(!stop.now){

        initial_weights <- c(W1, W2)

        curr_model = rgcca(A, weight.mat, tau = "optimal", verbose = FALSE,
            scheme = scheme)
        
        curr_g_score = as.matrix(A[[1]] %*% curr_model$a[[1]])
        curr_t_score = as.matrix(A[[2]] %*% curr_model$a[[2]])
        curr_p_score = as.matrix(A[[3]] %*% curr_model$a[[3]])
        
        curr_scores = cbind(curr_g_score, cbind(curr_t_score, curr_p_score))
        model_scores[[tx]] <- curr_scores

        #check signs before updating weights.
        curr_scores <- check_signs(curr_scores)
        
        #curr_cor = cor(curr_scores[[1]])
        if(use.partial.cor){
            curr_cor = pcor.shrink(curr_scores[[1]], verbose = FALSE)
        }else{
            curr_cor <- cor(curr_scores[[1]])
        }
        flag <- curr_scores[[2]]

        w1 = curr_cor[1, 2] / (1 - curr_cor[1, 2]^2)
        w2 = curr_cor[2, 3] / (1 - curr_cor[2, 3]^2)
        
        W1 = w1 / (w1 + w2)
        W2 = w2 / (w1 + w2)

        curr_weights <- c(W1, W2)
        if(verbose){print(c(W1, W2))}

        stopping.criteria <- check_stop(initial_weights, curr_weights, last_diff)
        last_diff <- stopping.criteria[[3]]
        stop.decision <- decide_stop(stopping.criteria, iter)
        stop.now <- stop.decision[[1]]
        iter = iter + 1
        
        weight.mat = 0.5 * matrix(c(0,0,0,W1,0,0,0,W2,0), 3, 3)
        weight.mat = weight.mat + t(weight.mat)
    }

    colnames(curr_scores[[1]]) <- c("Causal", "Mediator", "Outcome")
    curr_scores$reason <- stop.decision$reason
    return(curr_scores)
}
