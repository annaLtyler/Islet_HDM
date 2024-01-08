#xhist_prop and yhist_prop determine how much of the plotting space the x and
#y histograms take up; breaksx and breaksy determine the number of breaks for
#each histogram

plot_with_marginal_hist <- function(x,y, breaksx = 100, breaksy = 100, xhist_prop = 0.6,
    yhist_prop = 0.6, with.model = TRUE,
    xlim = NULL, ylim = NULL, col = "black", line.col = "#a6bddb", pch = 16, main, xlab, 
    ylab, report = c("lm", "cor.test"), cex = 1, plot.results = TRUE, add = FALSE, 
    plot.type = "p"){

	if(missing(main)){main = ""}
	if(missing(xlab)){xlab = deparse(substitute(x))}
	if(missing(ylab)){ylab = deparse(substitute(y))}

    report = report[1]

    layout.mat <- matrix(c(2,0,1,3), ncol = 2, byrow = TRUE)
    
    layout(layout.mat, widths = c(1,xhist_prop), heights = c(yhist_prop, 1))
    par(mar = c(4,4,0,0))
    if(with.model){
        result <- plot.with.model(x,y, xlab = xlab, ylab = ylab, xlim = xlim, 
            ylim = ylim, col = col, line.col = line.col, pch = pch, main = main, 
            report = report, cex = cex, plot.results = plot.results, add = add, 
            plot.type = plot.type, write.results = FALSE)
    }else{
        plot(x,y, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, cex = cex, 
            main = main)
    }
    plot.dim <- par("usr")
    par(mar = c(0,4,4,0))
    hist(x, axes = FALSE, breaks = breaksx, main = "", ylab = "")
    #axis(1)
    par(mar = c(4,0,0,4))
    bars <- hist(y, breaks = breaksy, plot = FALSE)
    barplot(bars$density, horiz = TRUE, axes = FALSE)
    #axis(2)

    if(with.model){
        if(report == "lm"){
            rlab <- "R2"
            }else{
            rlab <- "r"
        }

        new.title <- paste0(main, "\n", rlab, " = ", result[[1]], ", p = ", signif(result[[2]], 2))
    }else{
        new.title <- main
    }

        mtext(new.title, side = 3, outer = TRUE, line = -2.5)
    
    invisible(result)
}