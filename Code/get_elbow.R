#This function takes in a curve defined by x and y
#and returns its best guess for the elbow of the 
#curve. It returns the index of the point that is
#the farthest from the line defined by the first 
#point in x,y and the last point in x,y
#to plot where this elbow is on your curve, use
#plot.results = TRUE

get_elbow <- function(x,y, plot.results = FALSE){

    start.pt <- c(x[1],y[1])
    end.pt <- c(x[length(x)],y[length(y)])

    #for each point in the curve, find the distance to the
    #line defined by the endpoints.
    #find the slope and intercept of that line
    mb = intercept.slope(x0 = start.pt[1], x1 = end.pt[1], 
                         y0 = start.pt[2], y1 = end.pt[2])
    m = mb[1]
    b = mb[2]

    get.dist <- function(m, b, x0, y0){
        num <- abs(b + (m*x0) - y0)
        den <- sqrt((1 + m^2))
        distance <- num/den
        return(distance)
    }

    all.dist <- sapply(1:length(x), function(i) get.dist(m, b, x[i], y[i]))
    elbow.idx = which.max(all.dist)
    names(elbow.idx) <- "elbow_index"
    
    if(plot.results){
        plot(x,y)
        points(start.pt[1], start.pt[2], col = "red", pch = 16)
        points(end.pt[1], end.pt[2], col = "red", pch = 16)
        segments(start.pt[1], start.pt[2], end.pt[1], end.pt[2], col = "red")
        points(x[elbow.idx], y[elbow.idx], col = "red", pch = 16)    
    }
    
    return(elbow.idx)
}