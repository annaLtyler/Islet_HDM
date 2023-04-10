#This function takes in a list of matrices
#with rownames that have common IDs. The
#function find the common IDs across matrices
#and generates two new matrices that only 
#include common IDs and are in the same order.

get_common_rows <- function(matrix.list){
    common.ind <- Reduce("intersect", lapply(matrix.list, rownames))
    common.idx <- lapply(matrix.list, function(x) match(common.ind, rownames(x)))
    new.matrices <- lapply(1:length(matrix.list), function(x) matrix.list[[x]][common.idx[[x]],])
    return(new.matrices)
}