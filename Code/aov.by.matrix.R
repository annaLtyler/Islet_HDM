#This function assumes that each column of a matrix has
#its own factor. It creates a factor 
#based on the number of columns in the matrix.

aov.by.matrix <- function(matX, return.which = c("aov", "anova")){

    all.num <- as.vector(matX)
    fact.mat <- matrix(rep(1:ncol(matX), each = nrow(matX)), nrow = nrow(matX), ncol = ncol(matX))
    fact <- as.factor(as.vector(fact.mat))
    #boxplot(all.num~fact)
    if(return.which == "aov"){
        test <- aov(all.num~fact)
    }
    if(return.which == "anova"){
        test <- anova(aov(all.num~fact))
    }
    return(test)
}
