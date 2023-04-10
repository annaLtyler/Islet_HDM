#This function calculates the distance of each point
#in 2D space from the centroid of the distribution.
#It then ranks the points by distance to the center
dist_to_mid <- function(coord.mat, mid.fun = mean, plot.results = FALSE){

    mid_fun <- match.fun(mid.fun)
    mid.x <- mid_fun(coord.mat[,1], na.rm = TRUE)
    mid.y <- mid_fun(coord.mat[,2], na.rm = TRUE)

    pt.dist <- apply(coord.mat, 1, function(z) edist(z[1], mid.x, x[2], mid.y))
    dist.order <- order(pt.dist, decreasing = TRUE)

    if(plot.results){
        dist.col <- colors.from.values(pt.dist)
        plot(coord.mat, col = dist.col, pch = 16)
        points(mean.x, mean.y, pch = 16, col = "red")
    }
 
    return(coord.mat[dist.order,])

}
