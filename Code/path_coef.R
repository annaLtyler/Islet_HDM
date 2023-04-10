#In a high-dimensional mediation, you want the 
#correlation between the causal factor and the outcome
#to be equal to the product of the causal-mediator
#and mediator-outcome correlations.
#This function take in a correlation matrix in which
#the columns are ordered as causal, mediator, outcome
#and calculates the difference between the realized
#correlations and the ideal correlations.

path_coef <- function(model.scores){
    cor_mat <- pcor.shrink(model.scores, verbose = FALSE)
    path.coef <- cor_mat[1,2] * cor_mat[2,3]
    results <- c("path.coef" = path.coef, 
        "X-M" = cor_mat[1,2], "M-Y" = cor_mat[2,3],
        "X-Y" = cor_mat[1,3])
    return(results)
}
