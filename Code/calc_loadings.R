#This function calculates loadings on data
#after a high-dimensional mediation is performed
#with high_dim_med().
#The scores is the appropriate column from the
#scores matrix returned by high_dim_med(). Use
#drop=FALSE to make sure the scores are still a
#matrix. The data.mat is the data to compare to.
#The loadings are the correlation between the data
#and the scores. 
#For example, if you want to find the loadings onto
#transcripts that were used as mediators, use the
#Mediating column from the scores matrix, and the
#scaled expression matrix that went into the mediation.
#Make sure both have rownames that can be aligned
#by the function.

calc_loadings <- function(scores, data.mat){
    common.ind <- intersect(rownames(scores), rownames(data.mat))
    loadings = cor(data.mat[common.ind, ], scores, use = "pairwise.complete.obs")
    return(loadings)
}
