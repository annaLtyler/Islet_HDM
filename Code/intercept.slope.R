#find the intercept and slope of a line given x and y coordinates

intercept.slope <- function(x0, x1, y0, y1){
    m <- slope(x0, x1, y0, y1)
    b = y0 - (m*x0)
    result = c("slope" = m, "intercept" = b)
    return(result)
}
