#This function groups together values of one vector
#based on bins in another vector. 
#breaks is the same as in the function hist(). It
#can either be a positive integer or a vector of 
#values defining break points. This function draws
#barplots showing the mean of the binned values
#with standard errors. 

bin_by_hist <- function(x, y, breaks = 25, plot.result = TRUE,
    axes = TRUE, col = "gray", ylim = NULL, error.type = c("se", "sd"),
    xlab, ylab, main){

	if(missing(main)){main = ""}
	if(missing(xlab)){xlab = deparse(substitute(x))}
	if(missing(ylab)){ylab = deparse(substitute(y))}

    xhist <- hist(x, breaks = breaks, plot = FALSE)

    xmids <- xhist$mids
    xbins <- xhist$breaks
    
    yvals <- vector(mode = "list", length = length(xmids))
    for(i in 1:(length(xbins)-1)){
        bin.min <- xbins[i]
        bin.max <- xbins[(i+1)]
        x.idx <- intersect(which(x >= bin.min), which(x < bin.max))
        yvals[[i]] <- y[x.idx]
    }

    if(plot.result){
        ymeans <- sapply(yvals, function(x) mean(x, na.rm = TRUE))
        if(error.type == "se"){
            yse <- sapply(yvals, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
        }else{
            yse <- sapply(yvals, function(x) sd(x, na.rm = TRUE))
        }
        ymax <- max(c(ymeans, max(ymeans+yse, na.rm = TRUE)), na.rm = TRUE)
        if(is.null(ylim)){
            ylim = c(0, ymax)
        }
        a <- barplot(ymeans, ylim = ylim, col = col, axes = FALSE,
            xlab = xlab, ylab = ylab, main = main)
        axis(2, las = 2)
        segments(x0 = a, y0 = ymeans - yse, y1 = ymeans + yse, lwd = 3)
        if(axes){
            axis(1, at = a-0.5, labels = c(signif(xbins[1:length(xmids)], 2)))
        }
    result <- list("xmids" = xmids, "yvals" = yvals, "xcoord" = a)
    }else{
    result <- list("xmids" = xmids, "yvals" = yvals)
    }

    invisible(result)
}
