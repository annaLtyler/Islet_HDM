#This function gets peaks in a QTL table that overlap the input region.
#chr = 11; start.pos = 53.00257-2; end.pos = 53.00257+2
overlapping_peaks <- function(chr, start.pos, end.pos, peak.table){

    chr.locale <- which(peak.table[,"chr"] == chr)
    above.min <- which(peak.table[,"pos"] >= start.pos)
    below.max <- which(peak.table[,"pos"] <= end.pos)
    overlap.qtl <- Reduce("intersect", list(chr.locale, above.min, below.max))
    return(peak.table[overlap.qtl,])
}
