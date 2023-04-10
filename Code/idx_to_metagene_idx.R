#This function takes indices from a portrait matrix,
#(or a matrix of row/column numbers)
#in the orientation in which it is visualized, and
#translates them to metagene indices (which are in 
#a rotated matrix). R matrices start with 1,1 in the
#top left corner. Metagene matrices start with 1,1
#in the bottom left corner
idx_to_metagene_idx <- function(row.col = NULL, mat.idx = NULL, som.dim){

    if(!is.null(row.col)){
        mat.idx <- row_col_to_idx(row.col[,1], row.col[,2], som.dim)
    }

    new.mat <- matrix(0, nrow = som.dim, ncol = som.dim)
    new.mat[mat.idx] <- 1
    new.mat <- rotate.mat(new.mat)
    new.idx <- which(new.mat == 1)

    return(new.idx)

}
