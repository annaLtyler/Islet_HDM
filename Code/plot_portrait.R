#byrow indicates whether to build the matrix of values by row (TRUE)
#or by column (FALSE), and rotate is a numeric value 0, 1, 2, or 3
#indicating how many times the matrix should be rotated before displaying.

plot_portrait <- function(portraitV, global.min = NULL, global.max = NULL, main = "", 
  byrow = FALSE, rotate = 0, show.axes = FALSE, show.text = FALSE){
  
  portraitV[which(is.na(portraitV))] <- 0
  map.dim <- sqrt(length(portraitV))
  map.mat <- matrix(portraitV, nrow = map.dim, ncol = map.dim, byrow = byrow)
  if(rotate > 0){
    for(i in 1:rotate){
      map.mat <- rotate.mat(map.mat)
    }
  }
  if(!is.null(global.min)){
    imageWithText(map.mat, show.text = show.text, use.pheatmap.colors = TRUE, 
      global.color.scale = TRUE, global.min = global.min, global.max = global.max, 
      main = main)  
  }else{
    imageWithText(map.mat, show.text = show.text, use.pheatmap.colors = TRUE, main = main)
  }

  if(show.axes){
    axis(1, at = seq(0,ncol(map.mat), 5))
    axis(2, at = seq(0,nrow(map.mat), 5))
  }

  invisible(map.mat)
}
