#This function plots scan1 results from qtl2
#in a way that gives us more control over marking
#positions, like gene positions, etc.

plot_scan1 <- function(scan1.results, map, lodcol = 1, ylim = c(0, 8), 
    chr = c(1:19, "X"), main = "", type = "l", ylab = "LOD score",
    col = "black", mark.chr = NA, mark.pos = NA, mark.label = NA, 
    mark.font = 1, mark.col = "red", mark.lwd = 3, 
    mark.type = c("line", "arrow"), arrow.len = 0.25, arrow.max = 2, 
    arrow.min = 0, label.cex = 1){

    mark.type = mark.type[1]
    chr.x <- setup_genome_plot(map, ylim = ylim, chr = chr)
    plot.dim <- par("usr")
    plot.height <- plot.dim[4] - plot.dim[3]

    for(i in 1:length(map)){
        #if(i == 7){stop()}
        chr.pos <- map[[i]]
        marker.idx <- match(names(map[[i]]), rownames(scan1.results))
        x.pos <- scale.between.vals(chr.pos, chr.x[[i]][1], chr.x[[i]][2])
        points(x.pos, scan1.results[marker.idx,lodcol], type = type, col = col)

        if(names(chr.x)[i] %in% mark.chr){
            mark.idx <- which(mark.chr == names(chr.x)[i])
            mark.this <- mark.pos[mark.idx]
            nearest.pt <- get.nearest.pt(mark.this, map[[i]])

            if(mark.type == "line"){
                segments(x0 = x.pos[nearest.pt], y0 = (plot.dim[3]+(plot.height*0.03)), 
                    y1 = (plot.dim[4]-(plot.height*0.03)), col = mark.col, 
                    lwd = mark.lwd)
            }
            if(mark.type == "arrow"){
                arrows(x0 = x.pos[nearest.pt], y0 = arrow.max,
                y1 = arrow.min, len = arrow.len, col = mark.col, lwd = mark.lwd)
                if(!is.na(mark.label[mark.idx])){
                    text(x = x.pos[nearest.pt], y = (arrow.max+plot.height*0.03),
                    labels = mark.label[mark.idx], cex = label.cex, col = mark.col,
                    font = mark.font)
                }

            }
        }
    }

    axis(2)
    mtext(ylab, side = 2, line = 2.5)
    mtext(main, side = 3)
    mtext("Chromosome", side = 1, line = 0.5)
}