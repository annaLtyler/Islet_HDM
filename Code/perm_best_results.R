#This function takes results from get_perm_grid()
#and returns the best CCA penalty based on p 
#values and canonical correlations. If use.score is
#set to "cor.diff", the maximum difference between 
#the true and null canonical correlations is used.
#if use.score is set to "cor", the maximum canonical
#correlation is used, as long as it passes the p value
#threshold.
#the skew.toward parameters indicate whether to skew 
#toward harsher penalties (0), or more lenient penalties (1)
#if there are multiple best penalty pairs. The x and z
#allow skewing on both the x and the z axes independently.
#These parameters are only implemented if return.top.only
#is TRUE. Otherwise all best penalties are returned.

perm_best_results <- function(perm_grid, use.score = c("cor.diff", "cor"),
    pval.thresh = 0.05, return.top.only = FALSE, plot.results = TRUE, 
    row.text.shift = -0.25, col.text.shift = 0,
    skew.toward.x = 0, skew.toward.z = 0){

    use.score = use.score[1]
    pass.p <- which(perm_grid$p <= pval.thresh)
    if(use.score == "cor.diff"){
        cor.diff <- round(perm_grid$Cor - perm_grid$Perm.Cor, 2)
    }
    if(use.score == "cor"){
        cor.diff <- round(perm_grid$Cor, 2)
    }
    scale.factor <- 10000
    max.cor.diff <- floor(max(cor.diff[pass.p])*scale.factor)/scale.factor
    pass.cor <- which(cor.diff[pass.p] >= max.cor.diff)

    best.penalty.idx <- pass.p[pass.cor]

    if(length(best.penalty.idx) > 1 && return.top.only){ #pick the best penalty skewing toward harsher or more lenient depending on skew.toward.x and skew.toward.z
        penalty.row.col <- idx_to_row_col(best.penalty.idx, nrow(cor.diff))
        if(length(unique(penalty.row.col[,1])) > 1){
            best.row.idx <- get.nearest.pt(as.numeric(rownames(cor.diff)[penalty.row.col[,1]]), skew.toward.x)
        }else{
            best.row.idx <- 1
        }
        if(length(unique(penalty.row.col[,2])) > 1){
            best.col.idx <- get.nearest.pt(as.numeric(colnames(cor.diff)[penalty.row.col[,2]]), skew.toward.z)
        }else{
            best.col.idx <- 1
        }
        #pick an allowable best at random
        all.best <- sample(unique(c(best.row.idx, best.col.idx)), 1)
        best.penalty.idx <- row_col_to_idx(penalty.row.col[all.best, 1], penalty.row.col[all.best,2], nrow(cor.diff))
    }
    
    best.penalty.row.col <- idx_to_row_col(best.penalty.idx, nrow(cor.diff))
    
    best.penalty <- list("x" = as.numeric(rownames(cor.diff)[best.penalty.row.col[,1]]),
        "z" = as.numeric(colnames(cor.diff)[best.penalty.row.col[,2]]))

    star.nudge = 0.3

    if(plot.results){
        best.penalty.x <- best.penalty.row.col[1,1]
        best.penalty.z <- best.penalty.row.col[1,2]
        nx <- nrow(cor.diff)
        nz <- ncol(cor.diff)

        par(mfrow = c(2,2), mar = c(1,2,2,1))
        imageWithText(round(perm_grid$Cor, 2), col.scale = "blue", 
        main = "Correlations of Identified Components", 
        col.text.rotation = 0, row.text.shift = row.text.shift,
        col.text.shift = col.text.shift)
        points(x = best.penalty.z+star.nudge, 
        y = nx - best.penalty.x + 1 + star.nudge, 
        pch = "*", col = "red", cex = 1)
        
        imageWithText(round(perm_grid$Perm.Cor, 2), col.scale = "red", 
        main = "Correlations of Permuted Components", 
        col.text.rotation = 0, row.text.shift = row.text.shift,
        col.text.shift = col.text.shift)
        points(x = best.penalty.z+star.nudge, 
        y = nx - best.penalty.x + 1 + star.nudge, 
        pch = "*", col = "blue", cex = 1)

        imageWithText(round(cor.diff, 2), col.scale = "purple", 
        main = "Difference in Correlations", 
        col.text.rotation = 0, row.text.shift = row.text.shift,
        col.text.shift = col.text.shift)
        points(x = best.penalty.z+star.nudge, 
        y = nx - best.penalty.x + 1 + star.nudge, 
        pch = "*", col = "red", cex = 1)
        
        imageWithText(perm_grid$p, col.scale = "green", 
        main = "P values", col.text.rotation = 0, row.text.shift = row.text.shift,
        col.text.shift = col.text.shift)
        points(x = best.penalty.z+star.nudge, 
        y = nx - best.penalty.x + 1 + star.nudge, 
        pch = "*", col = "red", cex = 1)
    }

   return(best.penalty)

}