#This function returns points along a line 
#segment specified by x and y coordinates
#and a series of x or y coordinates to sample


get_line_pts <- function(x0, x1, y0, y1, sample.x = NULL, sample.y = NULL){
    
    if(length(sample.x) == 0 && length(sample.y) == 0){
        cat("Either sample.x or sample.y must be specified.")
        return(NULL)
    }

    if(length(sample.x) > 0 && length(sample.y) > 0){
        cat("Only on of sample.x or sample.y can be specified.")
        return(NULL)
    }

   mb <- intercept.slope(x0, x1, y0, y1)
   m = mb[1]
   b = mb[2]

    if(length(sample.x) > 0){
        y.pts <- sapply(sample.x, function(x) m*x + b)
        result <- cbind(sample.x, y.pts)
    }else{
        x.pts <- sapply(sample.y, function(y) (y-b)/m)
        result <- cbind(x.pts, sample.y)
    }
    colnames(result) <- c("x", "y")
    rownames(result) <- NULL
    return(result)

}