
#This function reassembles a matrix 
#so that when we use image to plot, 
#the matrix appears in the image in 
#the same orientation in which it
#is printed to the screen
#are in rows and the phenotypes in columns

rotate.mat <- function(mat){

	new.mat <- matrix(NA, nrow = dim(mat)[2], ncol = dim(mat)[1])
	#reverse each column i of the original
	#matrix and put it into row i in the
	#new matrix
	for(i in 1:dim(mat)[2]){
		new.mat[i,] <- rev(mat[,i])
		}
		rownames(new.mat) <- colnames(mat)
		colnames(new.mat) <- rev(rownames(mat))
		return(new.mat)
	}
