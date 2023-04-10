multilod.plot <- function(scan1.results, map, lod.thresh = 0, 
color.scheme = c("blue", "green", "purple", "red", "orange", "brown", "yellow", "gray"),
border.col = "darkgray", border.lwd = 3, row.names = NULL, row.name.shift = 0,
row.text.cex = 1, color.bar.cex = 1, color.bar.axis.lin = -2, 
color.fun = c("linear", "exponential"),
steepness = 1, mar = c(2, 6, 2, 0), chr.label.y = 0.5, chr.label.cex = 1, 
global.color.scale = FALSE,
global.min = lod.thresh, global.max = ceiling(max(scan1.results)), 
use.pheatmap.colors = FALSE){

    val.mat <- as.matrix(scan1.results)
    low.lod <- which(val.mat < lod.thresh)
    if(length(low.lod) > 0){
        val.mat[low.lod] <- NA
        has.vals <- which(apply(val.mat, 1, function(x) !all(is.na(x))))
        val.mat <- val.mat[has.vals,]
    }
    col.scale <- color.scheme[1]

    layout(matrix(c(1,2), ncol = 2), widths = c(1, 0.1))
    par(mar = mar)
    imageWithText(mat = t(val.mat), show.text = FALSE, col.names = NULL,
    col.scale = col.scale, color.fun = "exponential", exp.steepness = steepness,
    row.text.shift = row.name.shift, row.names = row.names, 
    row.text.cex = row.text.cex, global.color.scale = global.color.scale, 
    global.min = global.min, global.max = global.max, use.pheatmap.colors = use.pheatmap.colors)
    
    #label the chromosomes on the plot
    markers <- rownames(val.mat)
    all.markers <- lapply(map, names)
    
    for(i in 1:length(all.markers)){
        chr.locale <- which(markers %in% all.markers[[i]])
        chr.min <- min(chr.locale)
        chr.max <- max(chr.locale)
        chr.mid <- mean(c(chr.min, chr.max))
        draw.rectangle(chr.min, chr.max, 0.5, ncol(val.mat)+0.5, 
        border.col = border.col, lwd = border.lwd)
        text(x = chr.mid, y = chr.label.y, names(all.markers)[i], cex = chr.label.cex)
    }

    par(mar = c(2,2,2,2))
    bar.mat <- matrix(segment.region(min(val.mat, na.rm = TRUE), max(val.mat, na.rm = TRUE), 100), ncol = 1)
    imageWithTextColorbar(bar.mat, col.scale = col.scale, cex = color.bar.cex,
    axis.line = color.bar.axis.lin, color.fun = "exponential", exp.steepness = steepness,
    global.color.scale = TRUE, global.min = global.min, global.max = global.max, 
    use.pheatmap.colors = use.pheatmap.colors)

    invisible(val.mat)
}