hex2rgb <- function(hex.col, alpha = 1){
    rgb.vals <- col2rgb(hex.col)
    final.col <- rgb(rgb.vals[1]/256, rgb.vals[2]/256, rgb.vals[3]/256, alpha = alpha)
    return(final.col)
}
