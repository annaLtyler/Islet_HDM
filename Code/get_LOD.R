#This function finds the LOD score at a particular position

get_LOD <- function(scan1.obj, map, chr, pos){
    chr.locale <- which(names(map) == chr)
    nearest.marker <- get.nearest.pt(pos, map[[chr.locale]])
    marker.locale <- which(rownames(scan1.obj) == names(nearest.marker))
    marker.lod <- scan1.obj[marker.locale,]
    names(marker.lod) <- names(nearest.marker)
    return(marker.lod)
}