#This function finds the slope between adjacent points
#in x y vectors


instant.slope <- function(x, y){
	
	consec.mat <- consec_pairs(1:length(x))
	
	all.slopes <- apply(consec.mat, 1, function(v) slope(x[v[1]], x[v[2]], y[v[1]], y[v[2]]))
	
	return(all.slopes)

}