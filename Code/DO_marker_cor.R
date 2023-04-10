#chr = 11; start = 83; end = 85
DO_marker_cor <- function(genoprobs, map, chr, start, end){

    chr.idx <- which(names(genoprobs) == chr)
    chr.markers <- find_marker(map, chr = chr, interval = c(start, end))
    marker.pairs <- pair.matrix(1:length(chr.markers))
    
    marker.cor <- function(marker1, marker2){
        mcor <- cor(as.vector(genoprobs[[chr.idx]][,,marker1]), as.vector(genoprobs[[chr.idx]][,,marker2]))
        return(mcor)
    }

    cor.mat <- matrix(NA, nrow = length(chr.markers), ncol = length(chr.markers))
    rownames(cor.mat) <- colnames(cor.mat) <- chr.markers
    for(i in 1:nrow(marker.pairs)){
        cor.mat[marker.pairs[i,1], marker.pairs[i,2]] <- marker.cor(chr.markers[marker.pairs[i,1]], chr.markers[marker.pairs[i,2]])
    }
    #pheatmap(cor.mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)

    #marker.pair.cor <- apply(marker.pairs, 1, function(x) marker.cor(x[1], x[2]))

}