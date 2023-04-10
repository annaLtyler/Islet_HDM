#This function takes a matrix of numeric values and converts
#it into a 3D genotype array with the values scaled between
#0 and 1.

add_to_geno <- function(geno.obj, matX){

    #align the individuals in the two objects
    common.ind <- intersect(rownames(geno.obj), rownames(matX))
    common.geno.idx <- match(common.ind, rownames(geno.obj))
    common.mat.idx <- match(common.ind, rownames(matX))
    matched.geno <- geno.obj[common.geno.idx,,]
    matched.mat <- matX[common.mat.idx,,drop=FALSE]

  #scale to be between 0 and 1 for cape genotypes
  pos.mat <- apply(matched.mat, 2, function(x) x+abs(min(x)))
  #boxplot(pos.mat)
  scaled.mat <- apply(pos.mat, 2, function(x) x/max(x))
  #boxplot(scaled.mat)

  #create a 3D array of the scaled PCs
  new.geno <- array(NA, dim = c(nrow(matched.mat), ncol(matched.geno), ncol(matX)))
  dimnames(new.geno) <- list(rownames(matched.mat), colnames(matched.geno), colnames(scaled.mat))
  for(i in 1:ncol(scaled.mat)){
    one.mat <- matrix(scaled.mat[,i], nrow = nrow(scaled.mat), ncol = ncol(matched.geno))
    new.geno[,,i] <- one.mat
  }

  aug.geno <- abind(matched.geno, new.geno, along = 3)
  return(aug.geno)
}