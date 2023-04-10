#plot the partial correlation between two vectors
#in a set of vectors.
#mat is a matrix of values with features in 
#columns. 
#x.col is the number of the column to be plotted
#on the x axis and y.col is the number of the column
#to be plotted on the y axis.
#xlab and ylab are labels for the variables

plot.par.cor <- function(mat, x.col, y.col, xlab = colnames(mat)[x.col], 
ylab = colnames(mat)[y.col], main = ""){

    adj.pair <- adjust(mat[,c(x.col,y.col)], mat[,-c(x.col,y.col)])
    plot.with.model(adj.pair[,1], adj.pair[,2], xlab = xlab, ylab = ylab, 
    report = "cor.test", main = main)

}