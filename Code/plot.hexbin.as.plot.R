#This plots a hexbin plot using a regular plotting function so 
#we can get multiple plots on one page. I'm still not sure how
#do to it when plotting an actual hexbin object.

plot.hexbin.as.plot <- function(x, y, xlab, ylab, main, 
    min.cex = 1, max.cex = 3, n.bins = 10, count.scale.fun = NULL,
    legend.pos = "topright", round.legend = 10, use.pheatmap.colors = TRUE, 
    col.scale = "blue", grad.dir = "high", col.fun = c("linear", "exponential"), 
    exp.steepness = 1, light.dark = "f", custom.colors = NULL, legend.bg = par("bg"),
    with.model = FALSE, report = "lm", xlim = NULL, ylim = NULL){

    if(missing(main)){main = ""}
	if(missing(xlab)){xlab = deparse(substitute(x))}
	if(missing(ylab)){ylab = deparse(substitute(y))}
    col.fun <- col.fun[1]

    hb <- hexbin(x, y, xlab = xlab, ylab = ylab)
    #plot(hb)
    xy <- hcell2xy(hb)
    pt.cex = scale.between.vals(hb@count, target.min = min.cex, target.max = max.cex)
    if(!is.null(count.scale.fun)){
        count.scale.fun <- match.fun(count.scale.fun)
        cols <- colors.from.values(count.scale.fun(hb@count), 
        use.pheatmap.colors = use.pheatmap.colors,
        col.scale = col.scale, grad.dir = grad.dir, 
        color.fun = col.fun, exp.steepness = exp.steepness, light.dark = light.dark,
        custom.colors = custom.colors) 
    }else{
        cols <- colors.from.values(hb@count, use.pheatmap.colors = use.pheatmap.colors,
        col.scale = col.scale, grad.dir = grad.dir, color.fun = col.fun, light.dark = light.dark,
        exp.steepness = exp.steepness, custom.colors = custom.colors)    
    }
    
    plot.order = order(pt.cex, decreasing = FALSE)
    plot(lapply(xy, function(x) x[plot.order]), xlab = xlab, 
        ylab = ylab, pch = 16,
        main = main, col = cols[plot.order], 
        cex = pt.cex[plot.order], xlim = xlim, ylim = ylim)

    if(with.model){
        model <- lm(y~x)
        stats <- plot.with.model(x, y, report = report, plot.results = FALSE)
        abline(model)
        if(report == "lm"){
            stat.text <- paste0("R2 = ", signif(stats[1]), "; p = ", signif(stats[2]))
        }else{
            stat.text <- paste0("r = ", signif(stats[1]), "; p = ", signif(stats[2]))
        }
        mtext(stat.text, side = 3)
    }
    
    bins <- sort(unique(hb@count))
    if(is.null(round.legend)){
        round.legend <- max(bins)/3
    }
    binned.bins <- round(segment.region(min(bins), max(bins), n.bins, "ends"))
    rounded.bins <- unique(round(binned.bins/round.legend) * round.legend)
    if(length(rounded.bins) == 1){stop("I'm having trouble with the legend. Try reducing round.legend")}
    rounded.bins[which(rounded.bins == 0)] <- 1

    bin.cols <- colors.from.values(rounded.bins, use.pheatmap.colors = use.pheatmap.colors,
    col.scale = col.scale, grad.dir = grad.dir, color.fun = col.fun, exp.steepness = exp.steepness,
    light.dark = light.dark, custom.colors = custom.colors)

    bin.cex <- scale.between.vals(rounded.bins, target.min = min.cex, target.max = max.cex)
    legend(legend.pos, legend = rounded.bins, col = bin.cols, pch = 16, bg = legend.bg)

}