#This function is for plotting eQTLs and transcripts, etc.
#It supplies a relative plot position for a supplied 
#chromosome and position given a map. It allows longer
#chromosomes to appear longer in the plot.
#chr is a vector of chrommosomes
#pos is a vector of positions in Mb and is the same 
#length as chr
#map is the map for qtl2 mapping

get_relative_genomic_position <- function(chr, pos, map){
    
    #make relative positions for each SNP by chromosome
	 #get the max position for making relative positions of transcripts
	chr.max <- sapply(map, function(x) max(x))
	chr.sum <- sum(chr.max)

    chr.max.coord <- cumsum(chr.max)
    chr.min.coord <- chr.max.coord - chr.max

    #for each eQTL, figure out its relative position on the x axis
    rel.pos <- function(chr, pos){
        chr.locale <- which(names(map) == chr)
        rel.loc <- pos/chr.max[chr.locale]
        chr.len <- (chr.max.coord[chr.locale] - chr.min.coord[chr.locale])
        adjust.pos <- chr.len * rel.loc
        rel.x <-  chr.min.coord[chr.locale] + adjust.pos
        return(rel.x)
    }

    all.pos <- sapply(1:length(chr), function(x) rel.pos(chr[x], pos[x]))
    return(all.pos)
}