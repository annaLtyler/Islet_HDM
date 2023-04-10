#This function calculates the euclidean distance
#between two points specified by x and y coordinates.

edist <- function(x0, x1, y0, y1){
    the.dist <- sqrt((x0-x1)^2 + (y0 - y1)^2)
    return(the.dist)

}