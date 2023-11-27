#This function reorders a dendrogram so that it is
#maximall correlated with a numerical vector.
#The order corresponds to the rows of the original matrix.
reorder_dendrogram <- function(mat, numericV){

    d_row = dist(mat)
    hclust_row = hclust(d_row, method = "average")
    dend_row = as.dendrogram(hclust_row)
    dend_row = reorder(dend_row, numericV, agglo.FUN = median)
    order_dend_row = order.dendrogram(dend_row)
    #hclust_row = as.hclust(dend_row)
    return(order_dend_row)
}
