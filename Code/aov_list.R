#This function performs an ANOVA on a list.

aov_list <- function(listX){

    if(is.null(names(listX))){
        names(listX) <- paste0("Factor", 1:length(listX))
    }

    factor.names <- unlist(lapply(1:length(listX), function(x) rep(names(listX)[x], length(listX[[x]]))))
    vals <- unlist(listX)
    model <- aov(vals~as.factor(factor.names))
    return(model)
}